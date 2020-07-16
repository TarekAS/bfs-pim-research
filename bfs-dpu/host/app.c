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
  Row = 0,
  Col = 1,
  _2D = 2,
};

struct COO {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_nonzeros;
  uint32_t *row_idxs;
  uint32_t *col_idxs;
};

struct CSR {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_nonzeros;
  uint32_t row_offset;
  uint32_t *row_ptrs;
  uint32_t *col_idxs;
};

struct CSC {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_nonzeros;
  uint32_t *col_ptrs;
  uint32_t *row_idxs;
};

/**
 * @fn dpu_insert_mram_array_u32
 * @brief Inserts data into the MRAM of a DPU at the last used MRAM address.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol where to copy the pointer of the data.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_insert_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
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
 * @fn dpu_set_mram_array_u32
 * @brief Copy data to the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_set_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                      // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_array, (const uint8_t *)src, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_get_mram_array_u32
 * @brief Copy data from the MRAM of a DPU.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol of the pointer to MRAM destination.
 * @param dst the host buffer where the data is copied.
 * @param length the number of elements in the array.
 */
void dpu_get_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst, uint32_t length) {
  mram_addr_t p_array;
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                  // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_from_mram(dpu.dpu, (uint8_t *)dst, p_array, length * sizeof(uint32_t), 0)); // Copy data.
}

/**
 * @fn dpu_set_u32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param src the host buffer containing the data to copy
 */
void dpu_set_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t src) {
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &src, sizeof(uint32_t)));
}

/**
 * @fn dpu_get_u32
 * @brief Copy data from the Host memory buffer to one the DPU memories.
 * @param dpu_set the identifier of the DPU set
 * @param symbol_name the name of the DPU symbol where to copy the data
 * @param dst the host buffer where the data is copied
 */
void dpu_get_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *dst) {
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, dst, sizeof(uint32_t)));
}

// Parse CLI args and options.
void parse_args(int argc, char **argv, int *num_dpu, enum Algorithm *alg, enum Partition *prt, char **file) {
  bool is_prt_set = false;
  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "n:a:p:")) != -1)
    switch (c) {
    case 'n':
      *num_dpu = atoi(optarg);
      if (num_dpu == 0 || *num_dpu % 8 != 0) {
        PRINT_ERROR("Number of DPUs must be a multiple of 8.");
        exit(1);
      }
      break;
    case 'a':
      if (strcmp(optarg, "src") == 0) {
        *alg = SrcVtx;
        if (!is_prt_set)
          *prt = Row;
      } else if (strcmp(optarg, "dst") == 0) {
        *alg = DstVtx;
        if (!is_prt_set)
          *prt = Col;
      } else if (strcmp(optarg, "edge") == 0) {
        *alg = Edge;
        if (!is_prt_set)
          *prt = _2D;
      } else {
        PRINT_ERROR("Incorrect -a argument. Supported algorithms: src | dst | edge");
        exit(1);
      }
      break;
    case 'p':
      if (strcmp(optarg, "row") == 0)
        *prt = Row;
      else if (strcmp(optarg, "col") == 0)
        *prt = Col;
      else if (strcmp(optarg, "2d") == 0)
        *prt = _2D;
      else {
        PRINT_ERROR("Incorrect -p argument. Supported partitioning: row | col | 2d");
        exit(1);
      }
      is_prt_set = true;
      break;
    case '?':
    default:
      PRINT_ERROR("Bad args. Usage: -n <num_dpu> -a <src|dst|edge> -p <row|col|2d>");
      exit(1);
    }

  int numargs = argc - optind;
  if (numargs != 1) {
    if (numargs > 1)
      PRINT_ERROR("Too many arguments!");
    else
      PRINT_ERROR("Too few arguments! Please provide data file name (COO-formatted matrix).");
    exit(1);
  }
  *file = argv[optind];
}

// Load coo-formated file into memory.
// Pads nodes to multiple of lcm of 32 and n, and to a minimum of n * 32.
struct COO load_coo(char *file, uint32_t n) {

  if (access(file, F_OK) == -1) {
    PRINT_ERROR("Could not find file %s.", file);
    exit(1);
  }

  PRINT_INFO("Loading COO-formated graph from %s.", file);
  struct COO coo;

  // Initialize COO from file.
  uint32_t num_nodes = 0;
  uint32_t num_edges = 0;

  FILE *fp = fopen(file, "r");
  fscanf(fp, "%u", &num_nodes);
  fscanf(fp, "%u", &num_edges);
  coo.num_nonzeros = num_edges;
  coo.row_idxs = malloc(num_edges * sizeof(uint32_t));
  coo.col_idxs = malloc(num_edges * sizeof(uint32_t));

  // Pad number of nodes to a minimum of n * 32.
  uint32_t min = n * 32;
  if (num_nodes < min) {
    PRINT_WARNING("Number of nodes too low. Setting number of nodes to %d.", min);
    num_nodes = min;
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

  // Pad nodes to multiple of LCM.
  uint32_t padding = lcm - num_nodes % lcm;
  if (padding != 0) {
    PRINT_WARNING("Number of nodes must be multiple of %d. Padding with %d extra nodes.", lcm, padding);
    num_nodes += padding;
  }

  coo.num_rows = num_nodes;
  coo.num_cols = num_nodes;

  // Read nonzeros.
  PRINT_INFO("Reading COO-formated graph - %u nodes, %u edges.", num_nodes, num_edges);

  uint32_t row_idx, col_idx;
  fscanf(fp, "%u %u%*[^\n]\n", &row_idx, &col_idx); // Read first 2 integers of each line.

  uint32_t row_offset = row_idx; // Guarantee 0-indexed COO.
  coo.row_idxs[0] = row_idx - row_offset;
  coo.col_idxs[0] = col_idx - row_offset;

  for (uint32_t i = 1; i < num_edges; ++i) {
    fscanf(fp, "%u %u%*[^\n]\n", &row_idx, &col_idx);
    coo.row_idxs[i] = row_idx - row_offset;
    coo.col_idxs[i] = col_idx - row_offset;
  }
  fclose(fp);

  return coo;
}

// Partition COO matrix into n COO matrices by col, or by row, or both (2D). Assumes n is even.
struct COO *partition_coo(struct COO coo, int n, enum Partition prt) {

  struct COO *prts = malloc(n * sizeof(struct COO));

  // Initialize numNonzeros.
  for (int i = 0; i < n; ++i)
    prts[i].num_nonzeros = 0;

  uint32_t num_rows = coo.num_rows;
  uint32_t num_cols = coo.num_cols;
  uint32_t row_div = 1;
  uint32_t col_div = 1;

  // Determine num_rows, num_cols, and numNonZeros per partition.
  switch (prt) {
  case Row:
    num_rows /= n;
    for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
      uint32_t row_idx = coo.row_idxs[i];
      prts[row_idx / num_rows].num_nonzeros++;
    }
    break;

  case Col:
    num_cols /= n;
    for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
      uint32_t col_idx = coo.col_idxs[i];
      prts[col_idx / num_cols].num_nonzeros++;
    }
    break;

  case _2D:
    // Find the two nearest factors of n.
    row_div = (uint32_t)sqrt(n);
    while (n % row_div != 0)
      row_div--;
    col_div = n / row_div;

    num_rows /= row_div;
    num_cols /= col_div;

    for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
      uint32_t p_row = coo.row_idxs[i] / num_rows; // Partition row index.
      uint32_t p_col = coo.col_idxs[i] / num_cols; // Partition col index.
      uint32_t p = p_row * col_div + p_col;        // col-major index of coo.
      prts[p].num_nonzeros++;
    }

    break;
  }

  // Populate COO partitions.
  for (uint32_t i = 0; i < n; ++i) {
    prts[i].num_rows = num_rows;
    prts[i].num_cols = num_cols;
    prts[i].row_idxs = malloc(prts[i].num_nonzeros * sizeof(uint32_t));
    prts[i].col_idxs = malloc(prts[i].num_nonzeros * sizeof(uint32_t));
    prts[i].num_nonzeros = 0; // We'll reuse as index when appending data.
  }

  // Bin row and col pairs.
  for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    uint32_t col_idx = coo.col_idxs[i];

    uint32_t p;

    if (prt == Row)
      p = row_idx / num_rows;
    else if (prt == Col)
      p = col_idx / num_cols;
    else if (prt == _2D) {
      uint32_t p_row = coo.row_idxs[i] / num_rows;
      uint32_t p_col = coo.col_idxs[i] / num_cols;
      p = p_row * col_div + p_col;
    }

    uint32_t idx = prts[p].num_nonzeros;
    prts[p].col_idxs[idx] = col_idx;
    prts[p].row_idxs[idx] = row_idx;
    prts[p].num_nonzeros++;
  }

  return prts;
}

// Converts COO matrix to CSR format.
struct CSR coo_to_csr(struct COO coo) {

  struct CSR csr;

  // Initialize fields.
  csr.num_rows = coo.num_rows;
  csr.num_cols = coo.num_cols;
  csr.num_nonzeros = coo.num_nonzeros;
  csr.row_ptrs = calloc((csr.num_rows + 1), sizeof(uint32_t));
  csr.col_idxs = malloc(csr.num_nonzeros * sizeof(uint32_t));

  // Find smallest row_idx to use as offset. Usually it's the first row_idx, if data is sorted.
  uint32_t row_offset = coo.row_idxs[0];
  for (uint32_t i = 1; i < coo.num_nonzeros; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    if (row_idx < row_offset)
      row_offset = row_idx;
  }
  csr.row_offset = row_offset;

  // Histogram row_idxs.
  for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
    uint32_t row_idx = coo.row_idxs[i] - row_offset;
    csr.row_ptrs[row_idx]++;
  }

  // Prefix sum rowPtrs.
  uint32_t sum_before_next_row = 0;
  for (uint32_t row_idx = 0; row_idx < csr.num_rows; ++row_idx) {
    uint32_t sum_before_row = sum_before_next_row;
    sum_before_next_row += csr.row_ptrs[row_idx];
    csr.row_ptrs[row_idx] = sum_before_row;
  }
  csr.row_ptrs[csr.num_rows] = sum_before_next_row;

  // Bin the nonzeros.
  for (uint32_t i = 0; i < coo.num_nonzeros; ++i) {
    uint32_t row_idx = coo.row_idxs[i] - row_offset;
    uint32_t nnzIdx = csr.row_ptrs[row_idx]++;
    csr.col_idxs[nnzIdx] = coo.col_idxs[i];
  }

  // Restore rowPtrs.
  for (uint32_t row_idx = csr.num_rows - 1; row_idx > 0; --row_idx)
    csr.row_ptrs[row_idx] = csr.row_ptrs[row_idx - 1];
  csr.row_ptrs[0] = 0;

  return csr;
}

// Converts COO matrix to CSC format.
struct CSC coo_to_csc(struct COO coo) {

  // Transpose COO matrix.
  struct COO coo_trs = {
      .col_idxs = coo.row_idxs,
      .row_idxs = coo.col_idxs,
      .num_cols = coo.num_rows,
      .num_rows = coo.num_cols};

  // Convert to CSR, then CSC.
  struct CSR csr = coo_to_csr(coo_trs);
  struct CSC csc = {
      .col_ptrs = csr.row_ptrs,
      .row_idxs = csr.col_idxs,
      .num_cols = csr.num_cols,
      .num_rows = csr.num_rows};

  return csc;
}

void free_coo(struct COO coo) {
  free(coo.row_idxs);
  free(coo.col_idxs);
}

void free_csr(struct CSR csr) {
  free(csr.row_ptrs);
  free(csr.col_idxs);
}

void free_csc(struct CSC csc) {
  free(csc.col_ptrs);
  free(csc.row_idxs);
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

      dpu_get_mram_array_u32(dpu, "nextFrontier", nf_dpu, totalChunks);
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
    dpu_set_u32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_set_mram_array_u32(dpu, "nextFrontier", nextFrontier, totalChunks);
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
      dpu_get_u32(dpu, "numChunks", &numChunks);
      dpu_get_mram_array_u32(dpu, "nextFrontier",
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
    dpu_set_u32(set, "currentLevel", currentLevel);
    _DPU_FOREACH_I(set, dpu, i) {
      dpu_set_mram_array_u32(dpu, "currFrontier", nextFrontier, totalChunks);
    }
  }

  // Free resources.
  free(nextFrontier);
}

void print_node_levels(struct dpu_set_t set, struct dpu_set_t dpu, uint32_t numNodes) {
  // Get nodeLevels from each DPU.
  uint32_t *nodeLevels = malloc(numNodes * sizeof(uint32_t));
  uint32_t writeIdx = 0;

  int i = 0;
  _DPU_FOREACH_I(set, dpu, i) {
    uint32_t numNodes;
    dpu_get_u32(dpu, "numNodes", &numNodes);
    dpu_get_mram_array_u32(dpu, "nodeLevels", &nodeLevels[writeIdx], numNodes);
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

void start_src(struct dpu_set_t *set, struct dpu_set_t *dpu, int num_dpu, struct COO *coo_prts, enum Partition prt) {
  // Convert to CSR.
  struct CSR *csr = malloc(num_dpu * sizeof(struct CSR));
  for (int i = 0; i < num_dpu; ++i)
    csr[i] = coo_to_csr(coo_prts[i]);

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

  uint32_t i = 0;
  _DPU_FOREACH_I(*set, *dpu, i) {

    dpu_set_u32(*dpu, "dpu_idx", i);

    // Copy CSR partition.
    dpu_set_u32(*dpu, "num_nodes", csr[i].num_rows);
    dpu_set_u32(*dpu, "num_neighbors", csr[i].num_nonzeros);
    dpu_set_u32(*dpu, "node_offset", csr[i].row_offset);
    dpu_insert_mram_array_u32(*dpu, "node_ptrs", csr[i].row_ptrs, csr[i].num_rows + 1);
    dpu_insert_mram_array_u32(*dpu, "neighbors", csr[i].col_idxs, csr[i].num_nonzeros);

    PRINT_INFO("dpu_idx = %u, num_nodes = %u, num_neighbors = %u, node_offset = %u",
               i, csr[i].num_rows, csr[i].num_nonzeros, csr[i].row_offset);

    // dpu_insert_mram_array_u32(*dpu, "node_levels", 0, num_nodes);
    // dpu_insert_mram_array_u32(*dpu, "visited", 0, totalChunks);
    // dpu_insert_mram_array_u32(*dpu, "next_frontier", next_frontier, totalChunks);
    // dpu_insert_mram_array_u32(*dpu, "curr_frontier", 0, numChunks);
  }
}

int main(int argc, char **argv) {

  int num_dpu = 8;
  enum Algorithm alg = SrcVtx;
  enum Partition prt = Row;
  char *file = NULL;

  parse_args(argc, argv, &num_dpu, &alg, &prt, &file);

  char *bin_path;
  switch (alg) {
  case SrcVtx:
    bin_path = "bin/src-vtx";
    PRINT_INFO("Algorithm: Source-Vertex-based BFS.");
    break;
  case DstVtx:
    bin_path = "bin/dst-vtx";
    PRINT_INFO("Algorithm: Destination-Vertex-based BFS.");
    break;
  case Edge:
    bin_path = "bin/edge";
    PRINT_INFO("Algorithm: Edge-based BFS.");
    break;
  }

  switch (prt) {
  case Row:
    PRINT_INFO("1D Row partitioning (source-nodes).");
    break;
  case Col:
    PRINT_INFO("1D Column partitioning (destination-nodes/neighbors).");
    break;
  case _2D:
    PRINT_INFO("2D partitioning (both source-nodes and destination-nodes).");
    break;
  }

  struct dpu_set_t set, dpu;
  PRINT_INFO("Allocating %d DPUs.", num_dpu);
  DPU_ASSERT(dpu_alloc(num_dpu, NULL, &set));
  DPU_ASSERT(dpu_load(set, bin_path, NULL));

  struct COO coo = load_coo(file, num_dpu);                // Load coo-matrix from file.
  struct COO *coo_prts = partition_coo(coo, num_dpu, prt); // Partition COO by number of DPUs.
  free_coo(coo);

  switch (alg) {
  case SrcVtx:
    start_src(&set, &dpu, num_dpu, coo_prts, prt);
    break;

  case DstVtx: {
    // Convert to CSC.
    struct CSC *csc_prts = malloc(num_dpu * sizeof(struct CSR));
    for (int i = 0; i < num_dpu; ++i)
      csc_prts[i] = coo_to_csc(coo_prts[i]);
    // TODO: Populate MRAM (common alg)

    for (int i = 0; i < num_dpu; ++i)
      free_csc(csc_prts[i]);
    free(csc_prts);
  } break;

  case Edge:
    break;
  }

  for (int i = 0; i < num_dpu; ++i)
    free_coo(coo_prts[i]);
  free(coo_prts);

  return 0;

  print_node_levels(set, dpu, coo.num_rows);
  DPU_ASSERT(dpu_free(set));
  return 0;
}
