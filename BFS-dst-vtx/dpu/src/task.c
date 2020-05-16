#include <alloc.h>
#include <mram.h>
#include <stdint.h>
#include <stdio.h>
#include <seqread.h>
#include <defs.h>
#include <barrier.h>
#include <mutex.h>
#include <perfcounter.h>

#define PRINT_DEBUG(fmt, ...) printf("\033[0;34mDEBUG:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define ROUND_UP_TO_MULTIPLE_OF_8(x) ((((x) + 7) / 8) * 8)
#define ROUND_UP_TO_MULTIPLE_OF_32(x) ((((x)-1) / 32 + 1) * 32)

__host __mram_ptr void *p_used_mram_end = DPU_MRAM_HEAP_POINTER; // Points to the end of used MRAM addresses.

__host uint32_t numNodes;      // Number of nodes for this DPU.
__host uint32_t numNeighbors;  // Number of neighbors for this DPU.
__host uint32_t currentLevel;  // Current level in the BFS. Used to update the nodeLevel bitarray.
__host uint32_t totalChunks;   // Number of nodeChunks (i.e. length of nextFrontier).
__host uint32_t numChunks;     // Number of nodeChunks for this DPU (i.e. length of currFrontier). Should be divisible by NR_TASKLETS.

__host uint32_t nodeChunkFrom; // [nodeChunkFrom, nodeChunkTo[ is the range nodeChunks in nextFrontier for this DPU.
__host uint32_t nodeChunkTo;
__host uint32_t origin; // nodePtrs[0] of this DPU

__host __mram_ptr uint32_t *nodePtrs;        // DPU's share of nodePtrs.
__host __mram_ptr uint32_t *neighbors;       // DPU's share of neighbors.
__host __mram_ptr uint32_t *nextFrontier;    // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *currFrontier;    // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *visited;         // Nodes that are already visited.
__host __mram_ptr uint32_t *nodeLevels;      // The output of the BFS.

#ifndef NR_TASKLETS
#define NR_TASKLETS 16
#endif


BARRIER_INIT(nf_barrier, NR_TASKLETS);

int main() {

  // (void) perfcounter_config(COUNT_CYCLES, true);

  uint32_t chunksPerTasklet = numChunks / NR_TASKLETS;
  uint32_t idx = chunksPerTasklet * me();
  uint32_t lim = idx + chunksPerTasklet;

  // For each chunk of the currentFrontier belonging to this DPU. (Nested loop prevents false sharing)
  for (uint32_t c = idx; c < lim; ++c) {
    uint32_t cf = currFrontier[nodeChunkFrom + c];  // c: 0 -> numChunks

    nextFrontier[c] = 0;
    visited[c] |= cf;

    // For each (set) node in currFrontier, update its nodeLevel to the currentLevel.
    for (uint32_t b = 0; b < 32; ++b) {
      if(cf & (1 << (b % 32)))
        nodeLevels[c * 32 + b] = currentLevel;
    }
  }

  barrier_wait(&nf_barrier);

  // For each chunk of the nextFrontier belonging to this DPU. (Nested loop prevents false sharing)
  for (uint32_t c = idx; c < lim; ++c) {
    uint32_t unvis = ~visited[c];

    // For each unvisited node in the chunk.
    for (uint32_t b = 0; b < 32; ++b)
      if(unvis & (1 << (b % 32))) {
        uint32_t node = c * 32 + b;
        uint32_t offset = 1 << (node % 32);

        // Get nodePtrs of this node.
        uint32_t from = nodePtrs[node] - origin;
        uint32_t to;
        if (node < numNodes - 1)
          to = nodePtrs[node + 1] - origin;
        else
          to = numNeighbors;
        
        // For each neighbor.
        for (uint32_t n = from; n < to; ++n) {
          uint32_t neighbor = neighbors[n];
          uint32_t ncf = currFrontier[neighbor / 32]; // neighbor's currFrontier chunk.

          // If any neighbor is in currFrontier, add node to nextFrontier.
          if (ncf & (1 << (neighbor % 32))) {
            nextFrontier[c] |= offset;
            break;
          }
        }
      }
  }

  // perfcounter_t run_time = perfcounter_get();
  // printf("%lu ", run_time);
}
