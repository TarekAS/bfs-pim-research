#include <alloc.h>
#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <perfcounter.h>
#include <seqread.h>
#include <stdint.h>
#include <stdio.h>

#define PRINT_DEBUG(fmt, ...) printf("\033[0;34mDEBUG:\033[0m   " fmt "\n", ##__VA_ARGS__)

// Note: these are overriden by compiler flags.
#ifndef NR_TASKLETS
#define NR_TASKLETS 11
#endif
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif
#define BLOCK_INTS (BLOCK_SIZE / sizeof(uint32_t))
#ifndef BENCHMARK_CYCLES
#define BENCHMARK_CYCLES false
#endif

__host __mram_ptr void *p_used_mram_end = DPU_MRAM_HEAP_POINTER; // Points to the end of used MRAM addresses.

// CSC data.
__host __mram_ptr uint32_t *node_ptrs; // DPU's share of node_ptrs.
__host __mram_ptr uint32_t *edges;     // DPU's share of edges.

// BFS data.
__host uint32_t level;                     // Current level of the BFS.
__host uint32_t len_nf;                    // Length of next_frontier.
__host uint32_t nf_updated;                // DPU sets this to 1 if nf has been update in this level.
__host __mram_ptr uint32_t *visited;       // Nodes that are already visited.
__host __mram_ptr uint32_t *curr_frontier; // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *next_frontier; // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *node_levels;   // OUTPUT of the BFS.

// WRAM caches.
__dma_aligned uint32_t F_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t VIS_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t EDGE_CACHES[NR_TASKLETS][BLOCK_INTS];
__dma_aligned uint32_t NL_CACHES[NR_TASKLETS][32];

BARRIER_INIT(nf_barrier, NR_TASKLETS);

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
  uint32_t *edg = EDGE_CACHES[me()];
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

  // Loop over visited.
  for (uint32_t i = me() * BLOCK_INTS; i < len_nf; i += BLOCK_INTS * NR_TASKLETS) {
    mram_read(&visited[i], vis, BLOCK_SIZE);
    mram_read(&next_frontier[i], f, BLOCK_SIZE);

    for (uint32_t j = 0; j < BLOCK_INTS && i + j < len_nf; ++j) {

      uint32_t nonvis = ~vis[j];
      if (nonvis == 0)
        continue;

      // For each nonvisited node in the chunk.
      for (uint32_t b = 0; b < 32; ++b)
        if (nonvis & (1 << (b % 32))) {
          uint32_t node = (i + j) * 32 + b;
          uint32_t offset = 1 << (node % 32);

          // Get node_ptrs of this node.
          uint32_t from = node_ptrs[node];
          uint32_t to = node_ptrs[node + 1];

          // For each neighbor.
          for (uint32_t n = from; n < to; n += BLOCK_INTS) {

            // Handle &edges[from] not being 8-byte aligned.
            uint32_t k = 0;
            if (n == from && from % 2 != 0) {
              k = 1;
              n--;
            }
            mram_read(&edges[n], edg, BLOCK_SIZE);

            for (; k < BLOCK_INTS && n + k < to; ++k) {
              uint32_t neighbor = edg[k];
              uint32_t ncf = curr_frontier[neighbor / 32]; // neighbor's curr_frontier chunk.

              // If any neighbor is in curr_frontier, add node to next_frontier.
              if (ncf & (1 << (neighbor % 32))) {
                f[j] |= offset;
                nf_updated = 1;
                break;
              }
            }
          }
        }
    }
    mram_write(f, &next_frontier[i], BLOCK_SIZE);
  }
#if BENCHMARK_CYCLES
  cycles[me()] = perfcounter_get();
#endif
}
