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

// CSR data.
__host uint32_t num_nodes;             // Number of nodes for this DPU.
__host __mram_ptr uint32_t *node_ptrs; // DPU's share of node_ptrs.
__host __mram_ptr uint32_t *edges;     // DPU's share of edges.

// Chunks data.
__host uint32_t len_nf;  // Length of next_frontier.
__host uint32_t len_cf;  // Length of curr_frontier.
__host uint32_t cf_from; // [cf_from, cf_to[ ranges node chunks of this DPU.
__host uint32_t cf_to;

// BFS data.
__host uint32_t level;                     // Current level of the BFS.
__host __mram_ptr uint32_t *visited;       // Nodes that are already visited.
__host __mram_ptr uint32_t *curr_frontier; // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *next_frontier; // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *node_levels;   // OUTPUT of the BFS.

BARRIER_INIT(nf_barrier, NR_TASKLETS);
MUTEX_INIT(nf_mutex);

int main() {

  uint32_t nf_per_tasklet = len_nf / NR_TASKLETS;
  uint32_t idx = nf_per_tasklet * me();
  uint32_t lim = idx + nf_per_tasklet;

  // Loop over next_frontier.
  for (uint32_t c = idx; c < lim; ++c) {

    uint32_t f = next_frontier[c]; // Cache nf.
    next_frontier[c] = 0;          // Clear nf.
    visited[c] |= f;               // Update visited nodes.

    // Take DPU part of next_frontier as curr_frontier, and update node levels according to curr_frontier.
    if (c >= cf_from && c < cf_to) {
      uint32_t ridx = c - cf_from; // Relative index of curr_frontier in next_frontier.
      curr_frontier[ridx] = f;

      for (uint32_t b = 0; b < 32; ++b)
        if (f & 1 << b % 32)
          node_levels[ridx * 32 + b] = level;
    }
  }

  barrier_wait(&nf_barrier);

  // Loop over curr_frontier.
  for (uint32_t c = me(); c < len_cf; c += NR_TASKLETS) {

    uint32_t f = curr_frontier[c];

    // For each set node in the curr_frontier.
    for (uint32_t b = 0; b < 32; ++b)
      if (f & 1 << b % 32) {
        uint32_t node = c * 32 + b;

        // Get node_ptrs of this node.
        uint32_t from = node_ptrs[node];
        uint32_t to = node_ptrs[node + 1];

        // For each not visited neighbor of this node, add it to next_frontier.
        for (uint32_t n = from; n < to; ++n) {
          uint32_t neighbor = edges[n];
          uint32_t offset = 1 << neighbor % 32;

          if (!(visited[neighbor / 32] & offset)) {
            mutex_lock(nf_mutex);
            next_frontier[neighbor / 32] |= offset;
            mutex_unlock(nf_mutex);
          }
        }
      }
  }
}
