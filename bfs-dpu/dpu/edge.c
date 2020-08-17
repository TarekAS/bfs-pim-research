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

__host __mram_ptr void *p_used_mram_end = DPU_MRAM_HEAP_POINTER; // Points to the end of used MRAM addresses.

// COO data.
__host uint32_t num_edges;             // Length of nodes/dst_nodes.
__host uint32_t num_edges_tsk;         // Number of edges per tasklet.
__host __mram_ptr uint32_t *nodes;     // DPU's share of node idxs.
__host __mram_ptr uint32_t *neighbors; // DPU's share of neighbor idxs.

// Chunks data.
__host uint32_t len_nf;     // Length of next_frontier.
__host uint32_t len_nf_tsk; // Length of next_frontier per tasklet.

// BFS data.
__host uint32_t level;                     // Current level of the BFS.
__host __mram_ptr uint32_t *visited;       // Nodes that are already visited.
__host __mram_ptr uint32_t *curr_frontier; // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *next_frontier; // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *node_levels;   // OUTPUT of the BFS.

// Note: these are overriden by compiler flags.
#ifndef NR_TASKLETS
#define NR_TASKLETS 16
#endif

BARRIER_INIT(nf_barrier, NR_TASKLETS);
MUTEX_INIT(nf_mutex);

int main() {

  const uint32_t idx_nf = me() * len_nf_tsk;
  const uint32_t idx_edges = me() * num_edges_tsk;
  uint32_t lim_nf = idx_nf + len_nf_tsk;
  uint32_t lim_edges = idx_edges + num_edges_tsk;
  if (lim_edges > num_edges)
    lim_edges = num_edges;
  if (me() == NR_TASKLETS - 1)
    lim_nf = len_nf;

  // Loop over next_frontier.
  for (uint32_t c = idx_nf; c < lim_nf; ++c) {

    uint32_t f = next_frontier[c]; // Cache nf.
    visited[c] |= f;               // Update visited nodes.
    next_frontier[c] = 0;          // Clear nf.

    // Update node_levels according to the next_frontier.
    for (uint32_t b = 0; b < 32; ++b)
      if (f & (1 << (b % 32)))
        node_levels[c * 32 + b] = level;
  }

  barrier_wait(&nf_barrier);

  for (uint32_t edge = me(); edge < num_edges; edge += NR_TASKLETS) {
    uint32_t node = nodes[edge];
    if (curr_frontier[node / 32] & (1 << (node % 32))) {
      uint32_t neighbor = neighbors[edge];
      uint32_t offset = 1 << neighbor % 32;
      if (!(visited[neighbor / 32] & offset)) {
        mutex_lock(nf_mutex);
        next_frontier[neighbor / 32] |= offset;
        mutex_unlock(nf_mutex);
      }
    }
  }
}
