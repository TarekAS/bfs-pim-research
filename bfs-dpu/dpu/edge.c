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

__host uint32_t numNodes;      // Number of nodes for this DPU.
__host uint32_t numNeighbors;  // Number of neighbors for this DPU.
__host uint32_t currentLevel;  // Current level in the BFS. Used to update the nodeLevel bitarray.
__host uint32_t totalChunks;   // Number of nodeChunks (i.e. length of nextFrontier).
__host uint32_t numChunks;     // Number of nodeChunks for this DPU (i.e. length of currFrontier). Should be divisible by NR_TASKLETS.
__host uint32_t nodeChunkFrom; // [nodeChunkFrom, nodeChunkTo[ is the range nodeChunks in nextFrontier for this DPU.
__host uint32_t nodeChunkTo;
__host uint32_t origin;

__host __mram_ptr uint32_t *nodePtrs;     // DPU's share of nodePtrs.
__host __mram_ptr uint32_t *neighbors;    // DPU's share of neighbors.
__host __mram_ptr uint32_t *nextFrontier; // Nodes that are in the next frontier.
__host __mram_ptr uint32_t *currFrontier; // Nodes that are in the current frontier.
__host __mram_ptr uint32_t *visited;      // Nodes that are already visited.
__host __mram_ptr uint32_t *nodeLevels;   // The output of the BFS.

#ifndef NR_TASKLETS
#define NR_TASKLETS 16
#endif

BARRIER_INIT(nf_barrier, NR_TASKLETS);
MUTEX_INIT(nf_mutex);

int main() {

  uint32_t chunksPerTasklet = totalChunks / NR_TASKLETS;
  uint32_t idx = chunksPerTasklet * me();
  uint32_t lim = idx + chunksPerTasklet;

  // TODO: Implement edge-based BFS.

  for (uint32_t c = idx; c < lim; ++c) {
  }
  barrier_wait(&nf_barrier);

  for (uint32_t c = me(); c < numChunks; c += NR_TASKLETS) {
  }
}
