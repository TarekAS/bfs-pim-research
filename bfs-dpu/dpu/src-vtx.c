#include <alloc.h>
#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <perfcounter.h>
#include <seqread.h>
#include <stdint.h>
#include <stdio.h>

#define PRINT_DEBUG(fmt, ...) printf("\033[0;34mDEBUG:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define ROUND_UP_TO_MULTIPLE_OF_8(x) ((((x) + 7) / 8) * 8)
#define ROUND_UP_TO_MULTIPLE_OF_32(x) ((((x)-1) / 32 + 1) * 32)

// Note: this is overriden by compiler flags.
#ifndef NR_TASKLETS
#define NR_TASKLETS 16
#endif

__host __mram_ptr void *p_used_mram_end = DPU_MRAM_HEAP_POINTER; // Points to the end of used MRAM addresses.
__host uint32_t dpu_idx;                                         // DPU index.

// CSR data.
__host uint32_t num_nodes;             // Number of nodes for this DPU.
__host uint32_t num_neighbors;         // Number of neighbors for this DPU.
__host uint32_t node_offset;           // Offset of the first nodePtr relative to the global CSR.
__host __mram_ptr uint32_t *node_ptrs; // DPU's share of node_ptrs.
__host __mram_ptr uint32_t *neighbors; // DPU's share of neighbors.

// BFS metadata.
__host uint32_t level;        // Current level in the BFS.
__host uint32_t origin;       // Index of the first neighbor.
__host uint32_t num_chunks;   // Number of node chunks for this DPU (i.e. length of curr_frontier). Must be divisible by NR_TASKLETS.
__host uint32_t total_chunks; // Number of node chunks (i.e. length of next_frontier).
__host uint32_t chunk_from;   // [chunk_from, chunk_to[ ranges node chunks of this DPU.
__host uint32_t chunk_to;
__host __mram_ptr uint32_t *visited;       // Nodes that are already visited.
__host __mram_ptr uint32_t *curr_frontier; // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *next_frontier; // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *node_levels;   // OUTPUT of the BFS.

BARRIER_INIT(nf_barrier, NR_TASKLETS);
MUTEX_INIT(nf_mutex);

int main() {

  uint32_t chunksPerTasklet = total_chunks / NR_TASKLETS;
  uint32_t idx = chunksPerTasklet * me();
  uint32_t lim = idx + chunksPerTasklet;

  // For each chunk of the next_frontier.
  for (uint32_t c = idx; c < lim; ++c) {
    uint32_t f = next_frontier[c];

    next_frontier[c] = 0; // Clear chunk.
    visited[c] |= f;      // Update visited nodes.

    // If chunk belongs to this DPU.
    if (c >= chunk_from && c < chunk_to) {
      uint32_t ridx = c - chunk_from; // Relative index of next_frontier in the curr_frontier.
      curr_frontier[ridx] = f;        // Set curr_frontier.

      // For each (set) node in curr_frontier, update its nodeLevel to the level.
      for (uint32_t b = 0; b < 32; ++b)
        if (f & 1 << b % 32)
          node_levels[ridx * 32 + b] = level;
    }
  }
  barrier_wait(&nf_barrier);

  // For each chunk of the curr_frontier.
  for (uint32_t c = me(); c < num_chunks; c += NR_TASKLETS) {
    uint32_t f = curr_frontier[c];

    // For each set node in the curr_frontier.
    for (uint32_t b = 0; b < 32; ++b)
      if (f & 1 << b % 32) {
        uint32_t node = c * 32 + b;

        // Get node_ptrs of this node.
        uint32_t from = node_ptrs[node] - origin;
        uint32_t to;
        if (node < num_nodes - 1)
          to = node_ptrs[node + 1] - origin;
        else
          to = num_neighbors;

        // For each not visited neighbor of this node.
        for (uint32_t n = from; n < to; ++n) {
          uint32_t neighbor = neighbors[n];
          uint32_t offset = 1 << neighbor % 32;
          if (!(visited[neighbor / 32] & offset)) {

            // Add neighbor to next_frontier (critical section).
            mutex_lock(nf_mutex);
            next_frontier[neighbor / 32] |= offset;
            mutex_unlock(nf_mutex);
          }
        }
      }
  }
}
