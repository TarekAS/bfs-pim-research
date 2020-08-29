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

// Note: these are overriden by compiler flags.
#ifndef NR_TASKLETS
#define NR_TASKLETS 11
#endif
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif
#define BLOCK_INTS BLOCK_SIZE / sizeof(uint32_t)
#ifndef BENCHMARK_CYCLES
#define BENCHMARK_CYCLES false
#endif

__host __mram_ptr void *p_used_mram_end = DPU_MRAM_HEAP_POINTER; // Points to the end of used MRAM addresses.

// COO data.
__host uint32_t num_edges;             // Length of nodes/dst_nodes.
__host __mram_ptr uint32_t *nodes;     // DPU's share of node idxs.
__host __mram_ptr uint32_t *neighbors; // DPU's share of neighbor idxs.

// BFS data.
__host uint32_t level;                     // Current level of the BFS.
__host uint32_t len_nf;                    // Length of next_frontier.
__host __mram_ptr uint32_t *visited;       // Nodes that are already visited.
__host __mram_ptr uint32_t *curr_frontier; // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *next_frontier; // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *node_levels;   // OUTPUT of the BFS.

// WRAM caches.
__dma_aligned uint32_t F_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t VIS_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t NODES_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t NEIGHBORS_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t NL_CACHES[NR_TASKLETS][32];

BARRIER_INIT(nf_barrier, NR_TASKLETS);
MUTEX_INIT(nf_mutex);

#if BENCHMARK_CYCLES
__host uint64_t cycles[NR_TASKLETS];
#endif

int main() {
#if BENCHMARK_CYCLES
  if (me() == 0)
    (void)perfcounter_config(COUNT_CYCLES, true);
#endif

  uint32_t *f = F_CACHES[me()];
  uint32_t *vis = VIS_CACHES[me()];
  uint32_t *svtx = NODES_CACHES[me()];
  uint32_t *dvtx = NEIGHBORS_CACHES[me()];
  uint32_t *nl = NL_CACHES[me()];

  // Loop over next_frontier.
  for (uint32_t i = me() * BLOCK_INTS; i < len_nf; i += BLOCK_INTS * NR_TASKLETS) {
    mram_read(&visited[i], vis, BLOCK_SIZE);
    mram_read(&next_frontier[i], f, BLOCK_SIZE);

    for (uint32_t j = 0; j < BLOCK_INTS && i + j < len_nf; ++j) {
      uint32_t nf = f[j];
      if (nf == 0)
        continue;

      vis[j] |= nf; // Update visited nodes.
      f[j] = 0;     // Clear nf.

      // Update node levels.
      mram_read(&node_levels[(i + j) * 32], nl, 32 * sizeof(uint32_t));
      for (uint32_t b = 0; b < 32; ++b)
        if (nf & (1 << (b % 32)))
          nl[b] = level;
      mram_write(nl, &node_levels[(i + j) * 32], 32 * sizeof(uint32_t));
    }
    mram_write(vis, &visited[i], BLOCK_SIZE);
    mram_write(f, &next_frontier[i], BLOCK_SIZE);
  }

  barrier_wait(&nf_barrier);

  // Loop over edges.
  for (uint32_t i = me() * BLOCK_INTS; i < num_edges; i += BLOCK_INTS * NR_TASKLETS) {
    mram_read(&nodes[i], svtx, BLOCK_SIZE);
    mram_read(&neighbors[i], dvtx, BLOCK_SIZE);

    for (uint32_t j = 0; j < BLOCK_INTS && j + i < num_edges; ++j) {
      uint32_t node = svtx[j];
      if (curr_frontier[node / 32] & (1 << (node % 32))) {
        uint32_t neighbor = dvtx[j];
        uint32_t offset = 1 << neighbor % 32;
        if (!(visited[neighbor / 32] & offset)) {
          mutex_lock(nf_mutex);
          next_frontier[neighbor / 32] |= offset;
          mutex_unlock(nf_mutex);
        }
      }
    }
  }
#if BENCHMARK_CYCLES
  cycles[me()] = perfcounter_get();
#endif
}
