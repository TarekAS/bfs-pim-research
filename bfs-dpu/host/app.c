#include <assert.h>
#include <ctype.h>
#include <dpu.h>
#include <dpu_log.h>
#include <dpu_memory.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define _POSIX_C_SOURCE 2 // To use GNU's getopt.
#define PRINT_ERROR(fmt, ...) fprintf(stderr, "\033[0;31mERROR:\033[0m   " fmt "\n", ##__VA_ARGS__)
#define PRINT_WARNING(fmt, ...) fprintf(stderr, "\033[0;35mWARNING:\033[0m " fmt "\n", ##__VA_ARGS__)
#define PRINT_INFO(fmt, ...) fprintf(stderr, "\033[0;32mINFO:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define PRINT_STATUS(status) fprintf(stderr, "Status: %s\n", dpu_api_status_to_string(status))
#define PRINT_DEBUG(fmt, ...) printf("\033[0;34mDEBUG:\033[0m    " fmt "\n", ##__VA_ARGS__)

#define ROUND_UP_TO_MULTIPLE_OF_8(x) ((((x) + 7) / 8) * 8)
#define ROUND_UP_TO_MULTIPLE_OF_32(x) ((((x)-1) / 32 + 1) * 32)

enum Algorithm {
  SrcVtx = 0,
  DstVtx = 1,
  Edge = 2,
};

enum Partition {
  row = 0,
  col = 1,
  _2D = 2,
};

struct COOMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *rowIdxs;
  uint32_t *colIdxs;
};

struct CSRMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *rowPtrs;
  uint32_t *colIdxs;
};

struct CSCMatrix {
  uint32_t numRows;
  uint32_t numCols;
  uint32_t numNonzeros;
  uint32_t *colPtrs;
  uint32_t *rowIdxs;
};

// Parse CLI args and options.
void parse_args(int argc, char **argv, int *num_dpu, enum Algorithm *alg, enum Partition *prt, char **file) {
  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "n:a:p:")) != -1)
    switch (c) {
    case 'n':
      *num_dpu = atoi(optarg);
      if (num_dpu == 0 || *num_dpu % 8 != 0) {
        PRINT_ERROR("Number of DPUs must be a multiple of 8.\n");
        exit(1);
      }
      break;
    case 'a':
      if (strcmp(optarg, "src") == 0)
        *alg = SrcVtx;
      else if (strcmp(optarg, "dst") == 0)
        *alg = DstVtx;
      else if (strcmp(optarg, "edge") == 0)
        *alg = Edge;
      else {
        PRINT_ERROR("Incorrect -a argument. Supported algorithms: src | dst | edge\n");
        exit(1);
      }
      break;
    case 'p':
      if (strcmp(optarg, "row") == 0)
        *prt = row;
      else if (strcmp(optarg, "col") == 0)
        *prt = col;
      else if (strcmp(optarg, "2d") == 0)
        *prt = _2D;
      else {
        PRINT_ERROR("Incorrect -p argument. Supported partitioning: row | col | 2d\n");
        exit(1);
      }
      break;
    case '?':
    default:
      PRINT_ERROR("Bad args. Usage: -n <num_dpu> -a <src|dst|edge> -p <row|col|2d>\n");
      exit(1);
    }

  int numargs = argc - optind;
  if (numargs != 1) {
    if (numargs > 1)
      PRINT_ERROR("Too many arguments!\n");
    else
      PRINT_ERROR("Too few arguments! Please provide data file name (COO-formatted matrix).\n");
    exit(1);
  }
  *file = argv[optind];

  // Check if num_dpu is a perfect square in case of 2D partitioning.
  if (prt == _2D) {
    float fVar = sqrt((double)*num_dpu);
    int iVar = fVar;
    if (iVar != fVar) {
      PRINT_ERROR("Error: Number of DPUs must be a perfect square when choosing 2D partitioning.");
      exit(1);
    }
  }
}

// Load coo-formated file into memory.
// Pads rows/cols to multiple of lcm of 32 and n, and to a minimum of n * 32.
struct COOMatrix load_coo_matrix(char *file, uint32_t n) {

  if (access(file, F_OK) == -1) {
    PRINT_ERROR("Could not find file %s.", file);
    exit(1);
  }

  PRINT_INFO("Loading COO-formated graph from %s.", file);
  struct COOMatrix coo;

  // Initialize COO from file.
  FILE *fp = fopen(file, "r");
  fscanf(fp, "%u", &coo.numRows);
  fscanf(fp, "%u", &coo.numCols);
  fscanf(fp, "%u", &coo.numNonzeros);
  coo.rowIdxs = malloc(coo.numNonzeros * sizeof(uint32_t));
  coo.colIdxs = malloc(coo.numNonzeros * sizeof(uint32_t));

  if (coo.numRows != coo.numCols) {
    PRINT_ERROR("Error: number of rows and columns of COO-matrix are not equal. Exiting.");
    exit(1);
  }

  // Pad number rows/cols to a minimum of n * 32.
  uint32_t min = n * 32;
  if (coo.numRows < min) {
    PRINT_WARNING("Number of rows/cols too low. Setting number of rows/cols to %d.", min);
    coo.numRows = min;
    coo.numCols = min;
  }

  // Find Least Common Multiple of n and 32.
  uint32_t lcm = 0;
  uint32_t lar = (n > 32) ? n : 32;
  uint32_t small = (n > 32) ? 32 : n;

  for (int i = lar;; i += lar)
    if (i % small == 0) {
      lcm = i;
      break;
    }

  // Pad rows/cols to multiple of LCM.
  uint32_t padding = lcm - coo.numRows % lcm;

  if (padding != 0) {
    PRINT_WARNING("Number of rows/cols must be multiple of %d. Padding with %d extra rows/cols.", lcm, padding);
    coo.numRows += padding;
    coo.numCols += padding;
  }

  // Read nonzeros.
  PRINT_INFO("Reading COO-formated matrix - %u rows, %u columns, %u nonzeros.", coo.numRows, coo.numCols, coo.numNonzeros);

  uint32_t rowIdx, colIdx;
  fscanf(fp, "%u %u%*[^\n]\n", &rowIdx, &colIdx); // Read first 2 integers of each line.

  uint32_t nodeOffset = rowIdx; // Guarantee 0-indexed COO.
  coo.rowIdxs[0] = rowIdx - nodeOffset;
  coo.colIdxs[0] = colIdx - nodeOffset;

  for (uint32_t i = 1; i < coo.numNonzeros; ++i) {
    fscanf(fp, "%u %u%*[^\n]\n", &rowIdx, &colIdx);
    coo.rowIdxs[i] = rowIdx - nodeOffset;
    coo.colIdxs[i] = colIdx - nodeOffset;
  }
  fclose(fp);

  return coo;
}

// Partition COO matrix into n COO matrices by col, or by row, or both (2D).
struct COOMatrix *partition_coo(struct COOMatrix coo, int n, enum Partition prt) {

  struct COOMatrix *prts = malloc(n * sizeof(struct COOMatrix));

  // Initialize numNonzeros.
  for (int i = 0; i < n; ++i)
    prts[i].numNonzeros = 0;

  uint32_t rowDiv = 1;
  uint32_t colDiv = 1;

  // Determine numRows, numCols, and numNonZeros per partition.
  switch (prt) {
  case row:
    rowDiv = n; // n assumed to be even.
    for (int i = 0; i < coo.numNonzeros; ++i)
      prts[coo.rowIdxs[i] % n].numNonzeros++;
    break;

  case col:
    colDiv = n; // n assumed to be even.
    for (int i = 0; i < coo.numNonzeros; ++i)
      prts[coo.colIdxs[i] % n].numNonzeros++;
    break;

  case _2D:
    // Find the two nearest factors of n.
    rowDiv = (int)sqrt(n);
    while (n % rowDiv != 0)
      rowDiv--;
    colDiv = n / rowDiv;

    for (int i = 0; i < coo.numNonzeros; ++i) {
      uint32_t pRow = coo.rowIdxs[i] % rowDiv; // Partition row index.
      uint32_t pCol = coo.colIdxs[i] % colDiv; // Partition col index.
      uint32_t p = pRow * colDiv + pCol;       // col-major index of coo.
      prts[p].numNonzeros++;
    }

    break;
  }

  // Populate COO partitions.
  for (int i = 0; i < n; ++i) {
    prts[i].numRows = coo.numRows / rowDiv;
    prts[i].numCols = coo.numCols / colDiv;
    prts[i].rowIdxs = malloc(prts[i].numNonzeros * sizeof(uint32_t));
    prts[i].colIdxs = malloc(prts[i].numNonzeros * sizeof(uint32_t));
    prts[i].numNonzeros = 0; // We'll reuse as index to appending data.
  }

  // Bin row and col pairs.
  for (int i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    uint32_t colIdx = coo.colIdxs[i];

    uint32_t p;

    if (prt == row) {
      p = rowIdx % n;
    } else if (prt == col) {
      p = colIdx % n;
    } else if (prt == _2D) {
      uint32_t pRow = coo.rowIdxs[i] % (n / 2); // TODO: must be sqrt(n)
      uint32_t pCol = coo.colIdxs[i] % (n / 2); // TODO: must be sqrt(n)
      p = pRow * n + pCol;
    }

    uint32_t idx = prts[p].numNonzeros;
    prts[p].colIdxs[idx] = colIdx;
    prts[p].rowIdxs[idx] = rowIdx;
    prts[p].numNonzeros++;
  }

  return prts;
}

// Converts COO matrix to CSR format.
struct CSRMatrix coo_to_csr(struct COOMatrix coo) {

  struct CSRMatrix csr;

  // Initialize fields
  csr.numRows = coo.numRows;
  csr.numCols = coo.numCols;
  csr.numNonzeros = coo.numNonzeros;
  csr.rowPtrs = malloc((csr.numRows + 1) * sizeof(uint32_t));
  csr.colIdxs = malloc(csr.numNonzeros * sizeof(uint32_t));

  // Histogram rowIdxs
  memset(csr.rowPtrs, 0, (csr.numRows + 1) * sizeof(uint32_t));
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    csr.rowPtrs[rowIdx]++;
  }

  // Prefix sum rowPtrs
  uint32_t sumBeforeNextRow = 0;
  for (uint32_t rowIdx = 0; rowIdx < csr.numRows; ++rowIdx) {
    uint32_t sumBeforeRow = sumBeforeNextRow;
    sumBeforeNextRow += csr.rowPtrs[rowIdx];
    csr.rowPtrs[rowIdx] = sumBeforeRow;
  }
  csr.rowPtrs[csr.numRows] = sumBeforeNextRow;

  // Bin the nonzeros
  for (uint32_t i = 0; i < coo.numNonzeros; ++i) {
    uint32_t rowIdx = coo.rowIdxs[i];
    uint32_t nnzIdx = csr.rowPtrs[rowIdx]++;
    csr.colIdxs[nnzIdx] = coo.colIdxs[i];
  }

  // Restore rowPtrs
  for (uint32_t rowIdx = csr.numRows - 1; rowIdx > 0; --rowIdx) {
    csr.rowPtrs[rowIdx] = csr.rowPtrs[rowIdx - 1];
  }
  csr.rowPtrs[0] = 0;

  return csr;
}

// Converts COO matrix to CSC format.
struct CSCMatrix coo_to_csc(struct COOMatrix coo) {

  // Transpose COO matrix.
  struct COOMatrix cooTrs = coo;
  cooTrs.colIdxs = coo.rowIdxs;
  cooTrs.rowIdxs = coo.colIdxs;
  cooTrs.numCols = coo.numRows;
  cooTrs.numRows = coo.numCols;

  struct CSRMatrix csrMatrix = coo_to_csr(cooTrs);
  struct CSCMatrix csc = {
      .colPtrs = csrMatrix.rowPtrs,
      .rowIdxs = csrMatrix.colIdxs,
      .numCols = csrMatrix.numCols,
      .numRows = csrMatrix.numRows};

  return csc;
}

void free_coo_matrix(struct COOMatrix coo) {
  free(coo.rowIdxs);
  free(coo.colIdxs);
}

void free_csr_matrix(struct CSRMatrix csr) {
  free(csr.rowPtrs);
  free(csr.colIdxs);
}

void free_csc_matrix(struct CSCMatrix csc) {
  free(csc.colPtrs);
  free(csc.rowIdxs);
}

/**
 * @fn dpu_insert_to_mram_uint32_array
 * @brief Inserts data into the MRAM of a DPU at the last used MRAM address.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol where to copy the pointer of the data.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_insert_to_mram_uint32_array(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  // Get end of used MRAM pointer.
  mram_addr_t p_used_mram_end;
  DPU_ASSERT(dpu_copy_from(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Set the array pointer as the previous pointer.
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Copy the data to MRAM.
  size_t size = length * sizeof(uint32_t);
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_used_mram_end, (const uint8_t *)src, size, 0));

  // Increment end of used MRAM pointer.
  p_used_mram_end += size;
  DPU_ASSERT(dpu_copy_to(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));
}

/**
 * @fn dpu_copy_to_mram_uint32_array
 * @brief Copy data to the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_copy_to_mram_uint32_array(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                      // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_array, (const uint8_t *)src, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_copy_from_mram_uint32_array
 * @brief Copy data from the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param dst the host buffer where the data is copied.
 * @param length the number of elements in the array.
 */
void dpu_copy_from_mram_uint32_array(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                  // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_from_mram(dpu.dpu, (uint8_t *)dst, p_array, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_copy_to_uint32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param src the host buffer containing the data to copy
 */
void dpu_copy_to_uint32(struct dpu_set_t dpu, const char *symbol_name, uint32_t src) {
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &src, sizeof(uint32_t)));
}

/**
 * @fn dpu_copy_from_uint32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param dst the host buffer where the data is copied
 */
void dpu_copy_from_uint32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst) {
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, dst, sizeof(uint32_t)));
}

// Each DPU gets chunks of 32 nodes, distributed as evenly as possible.
uint32_t chunkize_nodes(uint32_t numRows, uint32_t num_dpu, uint32_t **dpu_node_chunks) {
  uint32_t totalChunks = (numRows + 31) / 32; // ceil.
  uint32_t minChunksPerDPU = totalChunks / num_dpu;
  uint32_t chunksRemainder = totalChunks % num_dpu;
  uint32_t *chunksPerDPU = malloc(num_dpu * sizeof(uint32_t));
  for (uint32_t i = 0; i < num_dpu; ++i)
    if (i < chunksRemainder)
      chunksPerDPU[i] = minChunksPerDPU + 1;
    else
      chunksPerDPU[i] = minChunksPerDPU;

  // Node idxs for each DPU.
  uint32_t *chunks = malloc((num_dpu + 1) * sizeof(uint32_t));
  chunks[0] = 0;
  for (uint32_t i = 1; i < num_dpu; ++i)
    chunks[i] = chunks[i - 1] + chunksPerDPU[i - 1] * 32;
  chunks[num_dpu] = numRows;

  *dpu_node_chunks = chunks;

  free(chunksPerDPU);
  return totalChunks;
}

void populate_mram(struct dpu_set_t *set, struct dpu_set_t *dpu, struct CSRMatrix csr, uint32_t *dpu_node_chunks, uint32_t totalChunks) {

  // Initialize nextFrontier.
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  nextFrontier[0] = 1; // Set root node.

  // Populate MRAM.
  PRINT_INFO("Populating MRAM.");
  dpu_copy_to_uint32(*set, "totalChunks", totalChunks);

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {
    uint32_t from = dpu_node_chunks[i];   // index of rowPtrs for this dpu.
    uint32_t to = dpu_node_chunks[i + 1]; // index of rowPtrs for the next dpu.
    uint32_t numNodes = to - from;
    uint32_t numNeighbors = csr.rowPtrs[to] - csr.rowPtrs[from];

    uint32_t nodeChunkFrom = from / 32;
    uint32_t nodeChunkTo = (to + 31) / 32;            // ceil.
    uint32_t numChunks = nodeChunkTo - nodeChunkFrom; // == chunksPerDPU[i]
    uint32_t *nodePtrs = &csr.rowPtrs[from];
    uint32_t origin = nodePtrs[0]; // index of first neighbor.
    uint32_t *neighbors = &csr.colIdxs[origin];

    PRINT_INFO("DPU %d: nodeChunks = [%u, %u[", i, nodeChunkFrom, nodeChunkTo);

    dpu_copy_to_uint32(*dpu, "numNodes", numNodes);
    dpu_copy_to_uint32(*dpu, "numNeighbors", numNeighbors);
    dpu_copy_to_uint32(*dpu, "numChunks", numChunks);
    dpu_copy_to_uint32(*dpu, "nodeChunkFrom", nodeChunkFrom);
    dpu_copy_to_uint32(*dpu, "nodeChunkTo", nodeChunkTo);
    dpu_copy_to_uint32(*dpu, "origin", origin);

    dpu_insert_to_mram_uint32_array(*dpu, "nodePtrs", nodePtrs, numNodes);
    dpu_insert_to_mram_uint32_array(*dpu, "neighbors", neighbors, numNeighbors);
    dpu_insert_to_mram_uint32_array(*dpu, "nodeLevels", 0, numNodes);

    dpu_insert_to_mram_uint32_array(*dpu, "visited", 0, totalChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "nextFrontier", nextFrontier, totalChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "currFrontier", 0, numChunks);
  }

  free(nextFrontier);
}

void populate_mram_csc(struct dpu_set_t *set, struct dpu_set_t *dpu, struct CSCMatrix csc, uint32_t *dpu_node_chunks, uint32_t totalChunks) {

  // Initialize nextFrontier.
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  nextFrontier[0] = 1; // Set root node.

  // Populate MRAM.
  PRINT_INFO("Populating MRAM.");
  dpu_copy_to_uint32(*set, "totalChunks", totalChunks);

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {
    uint32_t from = dpu_node_chunks[i];   // index of rowPtrs for this dpu.
    uint32_t to = dpu_node_chunks[i + 1]; // index of rowPtrs for the next dpu.
    uint32_t numNodes = to - from;
    uint32_t numNeighbors = csc.colPtrs[to] - csc.colPtrs[from];

    uint32_t nodeChunkFrom = from / 32;
    uint32_t nodeChunkTo = (to + 31) / 32;            // ceil.
    uint32_t numChunks = nodeChunkTo - nodeChunkFrom; // == chunksPerDPU[i]
    uint32_t *nodePtrs = &csc.colPtrs[from];
    uint32_t origin = nodePtrs[0]; // index of first neighbor.
    uint32_t *neighbors = &csc.rowIdxs[origin];

    PRINT_INFO("DPU %d: nodeChunks = [%u, %u[", i, nodeChunkFrom, nodeChunkTo);

    dpu_copy_to_uint32(*dpu, "numNodes", numNodes);
    dpu_copy_to_uint32(*dpu, "numNeighbors", numNeighbors);
    dpu_copy_to_uint32(*dpu, "numChunks", numChunks);
    dpu_copy_to_uint32(*dpu, "nodeChunkFrom", nodeChunkFrom);
    dpu_copy_to_uint32(*dpu, "nodeChunkTo", nodeChunkTo);
    dpu_copy_to_uint32(*dpu, "origin", origin);

    dpu_insert_to_mram_uint32_array(*dpu, "nodePtrs", nodePtrs, numNodes);
    dpu_insert_to_mram_uint32_array(*dpu, "neighbors", neighbors, numNeighbors);
    dpu_insert_to_mram_uint32_array(*dpu, "nodeLevels", 0, numNodes);

    dpu_insert_to_mram_uint32_array(*dpu, "visited", 0, numChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "nextFrontier", 0, numChunks);
    dpu_insert_to_mram_uint32_array(*dpu, "currFrontier", nextFrontier, totalChunks);
  }

  free(nextFrontier);
}

void bfs_src_vtx(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t totalChunks) {
  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  uint32_t *nf_dpu = calloc(totalChunks, sizeof(uint32_t));

  uint32_t currentLevel = 0;
  uint32_t done = true;
  int i = 0;

  while (true) {
    PRINT_INFO("Level %u", currentLevel);

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // Get union of nextFrontiers.
    _DPU_FOREACH_I(set, dpu, i) {

      dpu_copy_from_mram_uint32_array(dpu, "nextFrontier", nf_dpu, totalChunks);
      for (uint32_t c = 0; c < totalChunks; ++c) {
        nextFrontier[c] |= nf_dpu[c];
        if (nextFrontier[c] != 0)
          done = false;
      }

      // Get DPU logs.
      // PRINT_INFO("DPU %d:", i);
      // DPU_ASSERT(dpu_log_read(dpu, stdout));
    }

    if (done)
      break;

    done = true;
    ++currentLevel;

    // Update currentLevel and nextFrontier of DPUs.
    dpu_copy_to_uint32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_copy_to_mram_uint32_array(dpu, "nextFrontier", nextFrontier, totalChunks);
    }

    // Clear nextFrontier.
    for (uint32_t c = 0; c < totalChunks; ++c)
      nextFrontier[c] = 0;
  }

  free(nextFrontier);
  free(nf_dpu);
}

void bfs_dst_vtx(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t totalChunks) {

  uint32_t *nextFrontier = calloc(totalChunks, sizeof(uint32_t));
  uint32_t done = true;
  uint32_t currentLevel = 0;

  while (true) {
    PRINT_INFO("Level %u", currentLevel);

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

    // Concatenate all nextFrontiers.
    uint32_t writeIdx = 0;
    int i = 0;
    _DPU_FOREACH_I(set, dpu, i) {
      uint32_t numChunks;
      dpu_copy_from_uint32(dpu, "numChunks", &numChunks);
      dpu_copy_from_mram_uint32_array(dpu, "nextFrontier",
                                      &nextFrontier[writeIdx], numChunks);
      writeIdx += numChunks;

      // Get DPU logs.
      // PRINT_INFO("DPU %d", i);
      // DPU_ASSERT(dpu_log_read(dpu, stdout));
    }

    // Check if done.
    for (uint32_t c = 0; c < totalChunks; ++c)
      if (nextFrontier[c] != 0) {
        done = false;
        break;
      }
    if (done)
      break;

    done = true;
    ++currentLevel;

    // Update currentLevel and currentFrontier of DPUs.
    dpu_copy_to_uint32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_copy_to_mram_uint32_array(dpu, "currFrontier", nextFrontier, totalChunks);
    }
  }

  // Free resources.
  free(nextFrontier);
}

void bfs_edge(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t totalChunks) {
}

void print_node_levels(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t numNodes) {
  // Get nodeLevels from each DPU.
  uint32_t *nodeLevels = malloc(numNodes * sizeof(uint32_t));
  uint32_t writeIdx = 0;

  int i = 0;
  _DPU_FOREACH_I(set, dpu, i) {
    uint32_t numNodes;
    dpu_copy_from_uint32(dpu, "numNodes", &numNodes);
    dpu_copy_from_mram_uint32_array(dpu, "nodeLevels", &nodeLevels[writeIdx], numNodes);
    writeIdx += numNodes;
  }

  // Output.
  PRINT_INFO("Output:");
  for (uint32_t node = 0; node < numNodes; ++node) {
    uint32_t level = nodeLevels[node];
    if (node != 0 && level == 0) // Filters out "padded" rows.
      continue;
    printf("nodeLevels[%u]=%u\n", node, nodeLevels[node]);
  }

  free(nodeLevels);
}

int main(int argc, char **argv) {

  int num_dpu = 8;
  enum Algorithm alg = SrcVtx;
  enum Partition prt = row;
  char *file = NULL;

  parse_args(argc, argv, &num_dpu, &alg, &prt, &file);

  char *binPath;
  switch (alg) {
  case SrcVtx:
    binPath = "bin/src-vtx";
    break;
  case DstVtx:
    binPath = "bin/dst-vtx";
    break;
  case Edge:
    binPath = "bin/edge";
    break;
  }

  struct dpu_set_t set, dpu;
  PRINT_INFO("Allocating %d DPUs.", num_dpu);
  DPU_ASSERT(dpu_alloc(num_dpu, NULL, &set));
  DPU_ASSERT(dpu_load(set, binPath, NULL));

  struct COOMatrix coo = load_coo_matrix(file, num_dpu);         // Load coo-matrix from file.
  struct COOMatrix *coo_prts = partition_coo(coo, num_dpu, prt); // Partition COO by number of DPUs.
  free_coo_matrix(coo);

  if (alg == SrcVtx) {

    // Convert to CSR.
    struct CSRMatrix *csr_prts = malloc(num_dpu * sizeof(struct CSRMatrix));
    for (int i = 0; i < num_dpu; ++i)
      csr_prts[i] = coo_to_csr(coo_prts[i]);

    // TODO: Populate MRAM

    for (int i = 0; i < num_dpu; ++i)
      free_csr_matrix(csr_prts[i]); // Segfaulting. Are we initializing rowPtrs?
    free(csr_prts);

  } else if (alg == DstVtx) {

    // Convert to CSC.
    struct CSCMatrix *csc_prts = malloc(num_dpu * sizeof(struct CSRMatrix));
    for (int i = 0; i < num_dpu; ++i)
      csc_prts[i] = coo_to_csc(coo_prts[i]);
    // TODO: Populate MRAM (common alg)

    for (int i = 0; i < num_dpu; ++i)
      free_csc_matrix(csc_prts[i]);
    free(csc_prts);
  } else if (alg == Edge) {
    // TODO: Populate MRAM
  }

  // TODO: Launch specific algorithm.
  // TODO: More stuff

  for (int i = 0; i < num_dpu; ++i)
    free_coo_matrix(coo_prts[i]);
  free(coo_prts);

  print_node_levels(set, dpu, coo.numRows);
  DPU_ASSERT(dpu_free(set));
  return 0;
}
