#include <assert.h>
#include <ctype.h>
#include <dpu.h>
#include <dpu_log.h>
#include <dpu_memory.h>
#include <dpu_types.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define _POSIX_C_SOURCE 2 // To use GNU's getopt.
#define PRINT_ERROR(fmt, ...) fprintf(stderr, "\033[0;31mERROR:\033[0m   " fmt "\n", ##__VA_ARGS__)
#define PRINT_WARNING(fmt, ...) fprintf(stderr, "\033[0;35mWARN:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define PRINT_INFO(fmt, ...) fprintf(stderr, "\033[0;32mINFO:\033[0m    " fmt "\n", ##__VA_ARGS__)
#define PRINT_STATUS(status) fprintf(stderr, "Status: %s\n", dpu_api_status_to_string(status))
#define PRINT_DEBUG(fmt, ...) fprintf(stderr, "\033[0;34mDEBUG:\033[0m   " fmt "\n", ##__VA_ARGS__)
#define ROUND_UP_TO_MULTIPLE(x, y) ((((x)-1) / y + 1) * y)

// Note: these are overriden by compiler flags.
#ifndef NR_TASKLETS
#define NR_TASKLETS 11
#endif
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif
#ifndef BENCHMARK_CYCLES
#define BENCHMARK_CYCLES false
#endif
#ifndef BENCHMARK_TIME
#define BENCHMARK_TIME false
#endif

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
  uint32_t num_edges;
  uint32_t *row_idxs;
  uint32_t *col_idxs;
};

struct CSR {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_edges;
  uint32_t *row_ptrs;
  uint32_t *col_idxs;
};

struct CSC {
  uint32_t num_rows;
  uint32_t num_cols;
  uint32_t num_edges;
  uint32_t *col_ptrs;
  uint32_t *row_idxs;
};

FILE *out;
uint32_t num_dpu = 8;

struct dpu_symbol_t mram_heap_sym;
struct dpu_symbol_t level_sym;

mram_addr_t cf_addr;
mram_addr_t nf_addr;

struct dpu_set_t set;
struct dpu_set_t dpu;

// All timings are in seconds.
double dpu_compute_time = 0; // Total time spent in DPU computation.
double host_comm_time = 0;   // Total time spent by Host-DPU communication and data aggregation.
double pop_mram_time = 0;    // Time spent populating the MRAM (initial copy).
double fetch_res_time = 0;   // Time spent retrieving the results from MRAM (final copy).

/**
 * @fn dpu_insert_mram_array_u32
 * @brief Inserts data into the MRAM of a DPU at the last used MRAM address.
 * @param dpu_set the identifier of the DPU set.
 * @param symbol_name the name of the DPU symbol where to copy the pointer of the data.
 * @param src the host buffer containing the data to copy.
 * @param length the number of elements in the array.
 */
void dpu_insert_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {

  bool is_zero = src == 0;
  if (is_zero)
    src = calloc(length, sizeof(uint32_t));

  // Get end of used MRAM pointer.
  mram_addr_t p_used_mram_end;
  DPU_ASSERT(dpu_copy_from(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Set the array pointer as the previous pointer.
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Copy the data to MRAM.
  size_t size = length * sizeof(uint32_t);
  size += size % 8; // Guarantee address will be aligned on 8 bytes.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_used_mram_end, (const uint8_t *)src, size));

  // Increment end of used MRAM pointer.
  p_used_mram_end += size;
  DPU_ASSERT(dpu_copy_to(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));

  if (is_zero)
    free(src);
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
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));                   // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_array, (const uint8_t *)src, length * sizeof(uint32_t))); // Copy data.
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
  DPU_ASSERT(dpu_copy_from(dpu, symbol_name, 0, &p_array, sizeof(mram_addr_t)));               // Get address of array in MRAM.
  DPU_ASSERT(dpu_copy_from_mram(dpu.dpu, (uint8_t *)dst, p_array, length * sizeof(uint32_t))); // Copy data.
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

// Finds the two nearest factors of n.
void nearest_factors(uint32_t n, uint32_t *first, uint32_t *second) {
  uint32_t f = (uint32_t)sqrt(n);
  while (n % f != 0)
    f--;
  *first = f;
  *second = n / f;
}

// Parse CLI args and options.
void parse_args(int argc, char **argv, uint32_t *num_dpu, enum Algorithm *alg, enum Partition *prt, char **bin_path, char **file, char **out_file) {
  bool is_prt_set = false;
  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "n:a:p:o:")) != -1)
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
        PRINT_INFO("Algorithm: Source-Vertex-based BFS.");
        *bin_path = "bin/src-vtx-dma";
        *alg = SrcVtx;
        if (!is_prt_set)
          *prt = Row;
      } else if (strcmp(optarg, "dst") == 0) {
        PRINT_INFO("Algorithm: Destination-Vertex-based BFS.");
        *bin_path = "bin/dst-vtx-dma";
        *alg = DstVtx;
        if (!is_prt_set)
          *prt = Col;
      } else if (strcmp(optarg, "edge") == 0) {
        PRINT_INFO("Algorithm: Edge-based BFS.");
        *bin_path = "bin/edge-dma";
        *alg = Edge;
        if (!is_prt_set)
          *prt = _2D;
      } else {
        PRINT_ERROR("Incorrect -a argument. Supported algorithms: src | dst | edge");
        exit(1);
      }
      break;
    case 'p':
      if (strcmp(optarg, "row") == 0) {
        PRINT_INFO("Partitioning: 1D Row (source-nodes).");
        *prt = Row;
      } else if (strcmp(optarg, "col") == 0) {
        PRINT_INFO("Partitioning: 1D Column (destination-nodes/neighbors).");
        *prt = Col;
      } else if (strcmp(optarg, "2d") == 0) {
        PRINT_INFO("Partitioning: 2D (both source-nodes and destination-nodes).");
        *prt = _2D;
      } else {
        PRINT_ERROR("Incorrect -p argument. Supported partitioning: row | col | 2d");
        exit(1);
      }
      is_prt_set = true;
      break;
    case 'o':
      *out_file = optarg;
      break;
    case '?':
    default:
      PRINT_ERROR("Bad args. Usage: -n <num_dpu> -a <src|dst|edge> -p <row|col|2d> -o <output_file>");
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

  if (*out_file == NULL)
    *out_file = "/dev/null";
}

// Load coo-formated file into memory.
// Pads the number of nodes to guarantee divisibility by n and further divisibility by 32.
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
  int match = fscanf(fp, "%u %u", &num_nodes, &num_edges);
  if (match != 2) {
    PRINT_ERROR("Could not properly read COO file. First line must be of the form: NUM_NODES NUM_EDGES");
    exit(1);
  }

  coo.num_edges = num_edges;
  coo.row_idxs = malloc(num_edges * sizeof(uint32_t));
  coo.col_idxs = malloc(num_edges * sizeof(uint32_t));

  // Pad the number of nodes to guarantee divisibility by n and then by 32.
  uint32_t old = num_nodes;
  if (num_nodes % n != 0)
    num_nodes += n - num_nodes % n;

  uint32_t chunk_size = num_nodes / n;
  if (chunk_size % 32 != 0) {
    chunk_size += 32 - chunk_size % 32;
    num_nodes = chunk_size * n;
  }

  uint32_t padding = num_nodes - old;
  if (padding != 0)
    PRINT_WARNING("Padding number of nodes with %u extra nodes.", padding);

  coo.num_rows = num_nodes;
  coo.num_cols = num_nodes;

  // Read nonzeros.
  PRINT_INFO("Reading COO-formated graph - %u nodes, %u edges.", num_nodes, num_edges);

  uint32_t row_idx, col_idx;
  match = fscanf(fp, "%u %u%*[^\n]\n", &row_idx, &col_idx); // Read first 2 integers of file.
  if (match != 2) {
    PRINT_ERROR("Could not properly read line %u. Lines must be of the form: ROW_IDX COL_IDX", 1);
    exit(1);
  }

  uint32_t row_offset = row_idx; // Guarantee 0-indexed COO.
  coo.row_idxs[0] = row_idx - row_offset;
  coo.col_idxs[0] = col_idx - row_offset;

  for (uint32_t i = 1; i < num_edges; ++i) {
    match = fscanf(fp, "%u %u%*[^\n]\n", &row_idx, &col_idx);
    if (match != 2) {
      PRINT_ERROR("Could not properly read line %u. Lines must be of the form: ROW_IDX COL_IDX", i + 1);
      exit(1);
    }
    coo.row_idxs[i] = row_idx - row_offset;
    coo.col_idxs[i] = col_idx - row_offset;
  }
  fclose(fp);

  return coo;
}

// Partition COO matrix into n COO matrices by col, or by row, or both (2D). Assumes n is even.
struct COO *partition_coo(struct COO coo, uint32_t n, enum Partition prt) {

  PRINT_INFO("Partitioning COO-matrix into %u parts.", n);

  struct COO *prts = malloc(n * sizeof(struct COO));

  // Initialize num_edges.
  for (uint32_t i = 0; i < n; ++i)
    prts[i].num_edges = 0;

  uint32_t num_rows = coo.num_rows;
  uint32_t num_cols = coo.num_cols;
  uint32_t row_div = 1;
  uint32_t col_div = 1;
  bool offset_row = false;
  bool offset_col = false;

  // Determine num_rows, num_cols, and num_edges per partition.
  switch (prt) {
  case Row:
    offset_row = true;
    row_div = n;
    num_rows /= row_div;
    for (uint32_t i = 0; i < coo.num_edges; ++i) {
      uint32_t row_idx = coo.row_idxs[i];
      prts[row_idx / num_rows].num_edges++;
    }
    break;

  case Col:
    offset_col = true;
    col_div = n;
    num_cols /= col_div;
    for (uint32_t i = 0; i < coo.num_edges; ++i) {
      uint32_t col_idx = coo.col_idxs[i];
      prts[col_idx / num_cols].num_edges++;
    }
    break;

  case _2D:
    offset_row = true;
    offset_col = true;

    nearest_factors(n, &row_div, &col_div);
    num_rows /= row_div;
    num_cols /= col_div;

    for (uint32_t i = 0; i < coo.num_edges; ++i) {
      uint32_t p_row = coo.row_idxs[i] / num_rows; // Partition row index.
      uint32_t p_col = coo.col_idxs[i] / num_cols; // Partition col index.
      uint32_t p = p_row * col_div + p_col;        // col-major index of coo.
      prts[p].num_edges++;
    }
    break;
  }

  // Initialize COO partitions.
  for (uint32_t i = 0; i < n; ++i) {
    prts[i].num_rows = num_rows;
    prts[i].num_cols = num_cols;
    prts[i].row_idxs = malloc(prts[i].num_edges * sizeof(uint32_t));
    prts[i].col_idxs = malloc(prts[i].num_edges * sizeof(uint32_t));
    prts[i].num_edges = 0; // We'll re-increment as we append data.
  }

  // Bin row and col pairs.
  for (uint32_t i = 0; i < coo.num_edges; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    uint32_t col_idx = coo.col_idxs[i];

    uint32_t p = 0;

    if (prt == Row)
      p = row_idx / num_rows;
    else if (prt == Col)
      p = col_idx / num_cols;
    else if (prt == _2D) {
      uint32_t p_row = row_idx / num_rows;
      uint32_t p_col = col_idx / num_cols;
      p = p_row * col_div + p_col;
    }

    uint32_t idx = prts[p].num_edges;
    prts[p].row_idxs[idx] = row_idx;
    prts[p].col_idxs[idx] = col_idx;
    prts[p].num_edges++;
  }

  // Offset nodes.
  for (uint32_t p = 0; p < n; ++p) {
    uint32_t row_offset = offset_row ? p / col_div * num_rows : 0;
    uint32_t col_offset = offset_col ? p % col_div * num_cols : 0;

    for (uint32_t i = 0; i < prts[p].num_edges; ++i) {
      prts[p].row_idxs[i] -= row_offset;
      prts[p].col_idxs[i] -= col_offset;
    }
  }

  return prts;
}

// Converts COO matrix to CSR format.
struct CSR coo_to_csr(struct COO coo) {

  struct CSR csr;

  // Initialize fields.
  csr.num_rows = coo.num_rows;
  csr.num_cols = coo.num_cols;
  csr.num_edges = coo.num_edges;
  csr.row_ptrs = calloc((csr.num_rows + 1), sizeof(uint32_t));
  csr.col_idxs = malloc(csr.num_edges * sizeof(uint32_t));

  // Histogram row_idxs.
  for (uint32_t i = 0; i < coo.num_edges; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
    csr.row_ptrs[row_idx]++;
  }

  // Prefix sum row_ptrs.
  uint32_t sum_before_next_row = 0;
  for (uint32_t row_idx = 0; row_idx < csr.num_rows; ++row_idx) {
    uint32_t sum_before_row = sum_before_next_row;
    sum_before_next_row += csr.row_ptrs[row_idx];
    csr.row_ptrs[row_idx] = sum_before_row;
  }
  csr.row_ptrs[csr.num_rows] = sum_before_next_row;

  // Bin the nonzeros.
  for (uint32_t i = 0; i < coo.num_edges; ++i) {
    uint32_t row_idx = coo.row_idxs[i];
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
      .num_rows = coo.num_cols,
      .num_edges = coo.num_edges};

  // Convert to CSR, then CSC.
  struct CSR csr = coo_to_csr(coo_trs);
  struct CSC csc = {
      .col_ptrs = csr.row_ptrs,
      .row_idxs = csr.col_idxs,
      .num_cols = csr.num_rows,
      .num_rows = csr.num_cols,
      .num_edges = coo.num_edges};

  return csc;
}

// Frees COO matrix.
void free_coo(struct COO coo) {
  free(coo.row_idxs);
  free(coo.col_idxs);
}

// Frees CSR matrix.
void free_csr(struct CSR csr) {
  free(csr.row_ptrs);
  free(csr.col_idxs);
}

// Frees CSC matrix.
void free_csc(struct CSC csc) {
  free(csc.col_ptrs);
  free(csc.row_idxs);
}

// Prints the number of cycles of the worst performing DPU in the set.
void print_dpu_cycles(struct dpu_set_t set, struct dpu_set_t dpu) {
  uint64_t cycles[num_dpu][NR_TASKLETS];
  uint32_t i = 0;
  DPU_FOREACH(set, dpu, i) {
    DPU_ASSERT(dpu_prepare_xfer(dpu, &cycles[i]));
  }
  DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "cycles", 0, sizeof(uint64_t) * NR_TASKLETS, DPU_XFER_DEFAULT));

  // Get max cycles per DPU (among tasklets).
  uint64_t max_dpu_cycles[num_dpu];
  DPU_FOREACH(set, dpu, i) {
    uint32_t max = 0;
    for (uint32_t t = 0; t < NR_TASKLETS; t++) {
      uint64_t tasklet_cycles = cycles[i][t];
      if (tasklet_cycles > max)
        max = tasklet_cycles;
    }
    max_dpu_cycles[i] = max;
  }

  // Get avg and max DPU cycles per level (i.e. worst-performing DPU).
  uint64_t max_cycles_lvl = 0;
  for (uint32_t d = 0; d < num_dpu; ++d) {
    uint64_t max_dpu = max_dpu_cycles[d];
    if (max_dpu > max_cycles_lvl)
      max_cycles_lvl = max_dpu;
  }
  printf("%lu\n", max_cycles_lvl);
}

// Fetches and prints node levels from DPUs.
void print_node_levels(uint32_t total_nodes, uint32_t len_nl, uint32_t div) {
  fprintf(out, "node\tlevel\n");

#if BENCHMARK_TIME
  clock_t t = clock();
#endif

  uint32_t *node_levels = calloc(total_nodes, sizeof(uint32_t));
  uint32_t *nl_tmp = calloc(len_nl, sizeof(uint32_t));

  uint32_t i = 0;
  DPU_FOREACH(set, dpu, i) {
    dpu_get_mram_array_u32(dpu, "node_levels", nl_tmp, len_nl);
    for (uint32_t n = 0; n < len_nl; ++n) {
      uint32_t nreal = n + i / div * len_nl % total_nodes;
      if (node_levels[nreal] == 0 || nl_tmp[n] < node_levels[nreal])
        node_levels[nreal] = nl_tmp[n];
    }
  }

#if BENCHMARK_TIME
  t = clock() - t;
  fetch_res_time += ((double)t) / CLOCKS_PER_SEC;
#endif

  for (uint32_t node = 0; node < total_nodes; ++node) {
    uint32_t level = node_levels[node];
    if (node != 0 && level == 0) // Filters out "padded" rows.
      continue;
    fprintf(out, "%u\t%u\n", node, node_levels[node]);
  }

  free(nl_tmp);
  free(node_levels);
}

void bfs_vtx_row(uint32_t len_cf, uint32_t len_nf) {

  uint32_t size_nf = ROUND_UP_TO_MULTIPLE(len_nf * sizeof(uint32_t), 8);
  uint32_t size_cf = ROUND_UP_TO_MULTIPLE(len_cf * sizeof(uint32_t), 8);

  uint32_t *frontier = calloc(size_nf, 1);
  uint32_t *nf_tmp = calloc(size_nf * num_dpu, 1);
  uint32_t level = 0;
  bool done = true;

  while (true) {

#if BENCHMARK_TIME
    clock_t t = clock();
#endif

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

#if BENCHMARK_TIME
    t = clock() - t;
    dpu_compute_time += ((double)t) / CLOCKS_PER_SEC;
    t = clock();
#endif

    // Union next_frontiers.
    uint32_t i = 0;
    DPU_FOREACH(set, dpu, i) {
      DPU_ASSERT(dpu_prepare_xfer(dpu, &nf_tmp[i * len_nf]));
      // DPU_ASSERT(dpu_log_read(dpu, stderr));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_FROM_DPU, mram_heap_sym, nf_addr, size_nf, DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, i) {
      for (uint32_t c = 0; c < len_nf; ++c) {
        frontier[c] |= nf_tmp[i * len_nf + c];
        if (done && frontier[c] != 0)
          done = false;
      }
    }

#if BENCHMARK_CYCLES
    print_dpu_cycles(set, dpu);
#endif

    if (done)
      break;
    done = true;

    // Update level, next_frontier and current_frontier.
    ++level;
    DPU_ASSERT(dpu_copy_to_symbol(set, level_sym, 0, &level, sizeof(uint32_t)));
    DPU_ASSERT(dpu_prepare_xfer(set, frontier));
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_TO_DPU, mram_heap_sym, nf_addr, size_nf, DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, i) {
      DPU_ASSERT(dpu_prepare_xfer(dpu, &frontier[i * len_cf]));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_TO_DPU, mram_heap_sym, cf_addr, size_cf, DPU_XFER_DEFAULT));

    // Clear frontier.
    memset(frontier, 0, size_nf);

#if BENCHMARK_TIME
    t = clock() - t;
    host_comm_time += ((double)t) / CLOCKS_PER_SEC;
#endif
  }

  free(nf_tmp);
  free(frontier);
}

void bfs_vtx_col(uint32_t len_cf, uint32_t len_nf) {

  uint32_t size_nf = ROUND_UP_TO_MULTIPLE(len_nf * sizeof(uint32_t), 8);
  uint32_t size_cf = ROUND_UP_TO_MULTIPLE(len_cf * sizeof(uint32_t), 8);
  uint32_t *frontier = calloc(size_cf, 1);
  uint32_t level = 0;
  bool done = true;

  while (true) {

#if BENCHMARK_TIME
    clock_t t = clock();
#endif

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

#if BENCHMARK_TIME
    t = clock() - t;
    dpu_compute_time += ((double)t) / CLOCKS_PER_SEC;
    t = clock();
#endif

    // Concatenate all next_frontiers.
    uint32_t i = 0;
    DPU_FOREACH(set, dpu, i) {
      DPU_ASSERT(dpu_prepare_xfer(dpu, &frontier[i * len_nf]));
      // DPU_ASSERT(dpu_log_read(dpu, stderr));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_FROM_DPU, mram_heap_sym, nf_addr, size_nf, DPU_XFER_DEFAULT));

    // Check if done.
    for (uint32_t c = 0; c < len_cf; ++c)
      if (frontier[c] != 0) {
        done = false;
        break;
      }

#if BENCHMARK_CYCLES
    print_dpu_cycles(set, dpu);
#endif

    if (done)
      break;
    done = true;

    // Update level and curr_frontier of DPUs.
    ++level;
    DPU_ASSERT(dpu_copy_to_symbol(set, level_sym, 0, &level, sizeof(uint32_t)));
    DPU_ASSERT(dpu_prepare_xfer(set, frontier));
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_TO_DPU, mram_heap_sym, cf_addr, size_cf, DPU_XFER_DEFAULT));

#if BENCHMARK_TIME
    t = clock() - t;
    host_comm_time += ((double)t) / CLOCKS_PER_SEC;
#endif
  }

  free(frontier);
}

void bfs_vtx_2d(uint32_t len_frontier, uint32_t len_cf, uint32_t len_nf, uint32_t col_div) {

  uint32_t size_nf = ROUND_UP_TO_MULTIPLE(len_nf * sizeof(uint32_t), 8);
  uint32_t size_cf = ROUND_UP_TO_MULTIPLE(len_cf * sizeof(uint32_t), 8);
  uint32_t size_f = ROUND_UP_TO_MULTIPLE(len_frontier * sizeof(uint32_t), 8);
  uint32_t *frontier = calloc(size_f, 1);
  uint32_t *nf_tmp = calloc(size_nf * num_dpu, 1);
  uint32_t level = 0;
  bool done = true;

  while (true) {

#if BENCHMARK_TIME
    clock_t t = clock();
#endif

    // Launch DPUs.
    DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

#if BENCHMARK_TIME
    t = clock() - t;
    dpu_compute_time += ((double)t) / CLOCKS_PER_SEC;
    t = clock();
#endif

    // Concatenate by column and union by row the next_frontiers of each DPU.
    uint32_t i = 0;
    DPU_FOREACH(set, dpu, i) {
      DPU_ASSERT(dpu_prepare_xfer(set, &nf_tmp[i * len_nf]));
      // DPU_ASSERT(dpu_log_read(dpu, stderr));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_FROM_DPU, mram_heap_sym, nf_addr, size_nf, DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, i) {
      for (uint32_t c = 0; c < len_nf; ++c)
        frontier[c + i * len_nf % len_frontier] |= nf_tmp[i * len_nf + c];
    }

    // Check if done.
    for (uint32_t c = 0; c < len_frontier; ++c)
      if (frontier[c] != 0) {
        done = false;
        break;
      }

#if BENCHMARK_CYCLES
    print_dpu_cycles(set, dpu);
#endif

    if (done)
      break;
    done = true;

    // Update level, next_frontier and current_frontier.
    ++level;
    DPU_ASSERT(dpu_copy_to_symbol(set, level_sym, 0, &level, sizeof(uint32_t)));
    DPU_FOREACH(set, dpu, i) {
      DPU_ASSERT(dpu_prepare_xfer(dpu, &frontier[i * len_nf % len_frontier]));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_TO_DPU, mram_heap_sym, nf_addr, size_nf, DPU_XFER_DEFAULT));
    DPU_FOREACH(set, dpu, i) {
      DPU_ASSERT(dpu_prepare_xfer(dpu, &frontier[i / col_div * len_cf]));
    }
    DPU_ASSERT(dpu_push_xfer_symbol(set, DPU_XFER_TO_DPU, mram_heap_sym, cf_addr, size_cf, DPU_XFER_DEFAULT));

    // Clear frontier.
    memset(frontier, 0, size_f);

#if BENCHMARK_TIME
    t = clock() - t;
    host_comm_time += ((double)t) / CLOCKS_PER_SEC;
#endif
  }

  free(nf_tmp);
  free(frontier);
}

void start_src_vtx(struct COO *coo, int num_dpu, enum Partition prt) {

  // Convert COO partitions to CSR.
  struct CSR *csr = malloc(num_dpu * sizeof(struct CSR));
  for (int i = 0; i < num_dpu; ++i) {
    csr[i] = coo_to_csr(coo[i]);
    free_coo(coo[i]);
  }

  // Compute BFS metadata.
  uint32_t num_nodes = csr[0].num_rows;
  uint32_t num_neighbors = csr[0].num_cols;
  uint32_t len_cf = num_nodes / 32;
  uint32_t len_nf = num_neighbors / 32;
  uint32_t total_nodes, len_nl;
  uint32_t row_div = 1, col_div = 1;

  if (prt == Row) {
    row_div = num_dpu;
    total_nodes = num_neighbors;
    len_nl = num_nodes;
  } else if (prt == Col) {
    col_div = num_dpu;
    total_nodes = num_nodes;
    len_nl = num_nodes;
  } else {
    nearest_factors(num_dpu, &row_div, &col_div);
    total_nodes = num_nodes * num_dpu / col_div;
    len_nl = num_nodes;
  }

  uint32_t len_frontier = total_nodes / 32;
  uint32_t *frontier = calloc(len_frontier + BLOCK_SIZE, sizeof(uint32_t)); // +BLOCK_SIZE for safe margin of copy.
  frontier[0] = 1;                                                          // Set root node.

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

#if BENCHMARK_TIME
  clock_t t = clock();
#endif

  uint32_t i = 0;
  DPU_FOREACH(set, dpu, i) {

    // Copy BFS data.
    dpu_set_u32(dpu, "level", 0);
    dpu_set_u32(dpu, "len_nf", len_nf);
    dpu_set_u32(dpu, "len_cf", len_cf);

    // Add root node to cf of all DPUs of first row and to nf of all DPUs of first col.
    uint32_t *cf = i < col_div ? frontier : 0;
    uint32_t *nf = i % col_div == 0 ? frontier : 0;

    // Make sure arrays can be safely partitioned by NR_TASKLETS and BLOCK_SIZE.
    uint32_t lcf = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_cf, NR_TASKLETS), BLOCK_SIZE);
    uint32_t lnf = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_nf, NR_TASKLETS), BLOCK_SIZE);
    uint32_t lnl = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_nl, NR_TASKLETS), BLOCK_SIZE);

    dpu_insert_mram_array_u32(dpu, "visited", 0, lnf);
    dpu_insert_mram_array_u32(dpu, "next_frontier", nf, lnf);
    dpu_insert_mram_array_u32(dpu, "curr_frontier", cf, lcf);
    dpu_insert_mram_array_u32(dpu, "node_levels", 0, lnl);

    // Copy CSR data. Variable sized buffers must be copied last.
    dpu_insert_mram_array_u32(dpu, "node_ptrs", csr[i].row_ptrs, num_nodes + 1);
    dpu_insert_mram_array_u32(dpu, "edges", csr[i].col_idxs, csr[i].num_edges);

    // Cache some MRAM addresses (address must be the same for all DPUs).
    DPU_ASSERT(dpu_copy_from(dpu, "next_frontier", 0, &nf_addr, sizeof(mram_addr_t)));
    DPU_ASSERT(dpu_copy_from(dpu, "curr_frontier", 0, &cf_addr, sizeof(mram_addr_t)));
  }

#if BENCHMARK_TIME
  t = clock() - t;
  pop_mram_time += ((double)t) / CLOCKS_PER_SEC;
#endif

  // Free resources.
  free(frontier);
  for (int i = 0; i < num_dpu; ++i)
    free_csr(csr[i]);

  // Start BFS algorithm.
  PRINT_INFO("Starting BFS algorithm.");
  if (prt == Row)
    bfs_vtx_row(len_cf, len_nf);
  else if (prt == Col)
    bfs_vtx_col(len_cf, len_nf);
  else
    bfs_vtx_2d(len_frontier, len_cf, len_nf, col_div);

  // Print node levels.
  print_node_levels(total_nodes, len_nl, col_div);
}

void start_dst_vtx(struct COO *coo, int num_dpu, enum Partition prt) {

  // Convert COO partitions to CSC.
  struct CSC *csc = malloc(num_dpu * sizeof(struct CSC));
  for (int i = 0; i < num_dpu; ++i) {
    csc[i] = coo_to_csc(coo[i]);
    free_coo(coo[i]);
  }

  // Compute BFS metadata.
  uint32_t num_nodes = csc[0].num_rows;
  uint32_t num_neighbors = csc[0].num_cols;
  uint32_t len_cf = num_nodes / 32;
  uint32_t len_nf = num_neighbors / 32;
  uint32_t total_nodes, len_nl;
  uint32_t row_div = 1, col_div = 1;

  if (prt == Row) {
    row_div = num_dpu;
    total_nodes = num_neighbors;
    len_nl = num_neighbors;
  } else if (prt == Col) {
    col_div = num_dpu;
    total_nodes = num_nodes;
    len_nl = num_neighbors;
  } else {
    nearest_factors(num_dpu, &row_div, &col_div);
    total_nodes = num_nodes * num_dpu / col_div;
    len_nl = num_neighbors;
  }

  uint32_t len_frontier = total_nodes / 32;
  uint32_t *frontier = calloc(len_frontier + BLOCK_SIZE, sizeof(uint32_t)); // +BLOCK_SIZE for safe margin of copy.
  frontier[0] = 1;                                                          // Set root node.

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

#if BENCHMARK_TIME
  clock_t t = clock();
#endif

  uint32_t i = 0;
  DPU_FOREACH(set, dpu, i) {

    // Copy BFS data.
    dpu_set_u32(dpu, "level", 0);
    dpu_set_u32(dpu, "len_nf", len_nf);

    // Make sure arrays can be safely partitioned by NR_TASKLETS and BLOCK_SIZE.
    uint32_t lcf = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_cf, NR_TASKLETS), BLOCK_SIZE);
    uint32_t lnf = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_nf, NR_TASKLETS), BLOCK_SIZE);
    uint32_t lnl = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_nl, NR_TASKLETS), BLOCK_SIZE);

    // Add root node to cf of all DPUs of first row and to nf of all DPUs of first col.
    uint32_t *cf = i < col_div ? frontier : 0;
    uint32_t *nf = i % col_div == 0 ? frontier : 0;

    dpu_insert_mram_array_u32(dpu, "visited", 0, lnf);
    dpu_insert_mram_array_u32(dpu, "next_frontier", nf, lnf);
    dpu_insert_mram_array_u32(dpu, "curr_frontier", cf, lcf);
    dpu_insert_mram_array_u32(dpu, "node_levels", 0, lnl);

    // Copy CSR data. Variable sized buffers must be copied last.
    dpu_insert_mram_array_u32(dpu, "node_ptrs", csc[i].col_ptrs, num_neighbors + 1);
    dpu_insert_mram_array_u32(dpu, "edges", csc[i].row_idxs, csc[i].num_edges);

    // Cache some MRAM addresses (address must be the same for all DPUs).
    DPU_ASSERT(dpu_copy_from(dpu, "next_frontier", 0, &nf_addr, sizeof(mram_addr_t)));
    DPU_ASSERT(dpu_copy_from(dpu, "curr_frontier", 0, &cf_addr, sizeof(mram_addr_t)));
  }

#if BENCHMARK_TIME
  t = clock() - t;
  pop_mram_time += ((double)t) / CLOCKS_PER_SEC;
#endif

  // Free resources.
  free(frontier);
  for (int i = 0; i < num_dpu; ++i)
    free_csc(csc[i]);

  // Start BFS algorithm.
  PRINT_INFO("Starting BFS algorithm.");
  if (prt == Row)
    bfs_vtx_row(len_cf, len_nf);
  else if (prt == Col)
    bfs_vtx_col(len_cf, len_nf);
  else
    bfs_vtx_2d(len_frontier, len_cf, len_nf, col_div);

  // Print node levels.
  print_node_levels(total_nodes, len_nl, 1);
}

void start_edge(struct COO *coo, int num_dpu, enum Partition prt) {

  // Compute BFS metadata.
  uint32_t num_nodes = coo[0].num_rows;
  uint32_t num_neighbors = coo[0].num_cols;
  uint32_t len_cf = num_nodes / 32;
  uint32_t len_nf = num_neighbors / 32;
  uint32_t total_nodes, len_nl;
  uint32_t row_div = 1, col_div = 1;

  if (prt == Row) {
    row_div = num_dpu;
    total_nodes = num_neighbors;
    len_nl = num_neighbors;
  } else if (prt == Col) {
    col_div = num_dpu;
    total_nodes = num_nodes;
    len_nl = num_neighbors;
  } else {
    nearest_factors(num_dpu, &row_div, &col_div);
    total_nodes = num_nodes * num_dpu / col_div;
    len_nl = num_neighbors;
  }

  uint32_t len_frontier = total_nodes / 32;
  uint32_t *frontier = calloc(len_frontier + BLOCK_SIZE, sizeof(uint32_t)); // +BLOCK_SIZE for safe margin of copy.
  frontier[0] = 1;                                                          // Set root node.

  // Copy data to MRAM.
  PRINT_INFO("Populating MRAM.");

#if BENCHMARK_TIME
  clock_t t = clock();
#endif

  uint32_t i = 0;
  DPU_FOREACH(set, dpu, i) {

    // Copy BFS data.
    dpu_set_u32(dpu, "level", 0);
    dpu_set_u32(dpu, "len_nf", len_nf);

    // Add root node to cf of all DPUs of first row and to nf of all DPUs of first col.
    uint32_t *cf = i < col_div ? frontier : 0;
    uint32_t *nf = i % col_div == 0 ? frontier : 0;

    // Make sure arrays can be safely partitioned by NR_TASKLETS and BLOCK_SIZE.
    uint32_t lcf = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_cf, NR_TASKLETS), BLOCK_SIZE);
    uint32_t lnf = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_nf, NR_TASKLETS), BLOCK_SIZE);
    uint32_t lnl = ROUND_UP_TO_MULTIPLE(ROUND_UP_TO_MULTIPLE(len_nl, NR_TASKLETS), BLOCK_SIZE);

    dpu_insert_mram_array_u32(dpu, "visited", 0, lnf);
    dpu_insert_mram_array_u32(dpu, "next_frontier", nf, lnf);
    dpu_insert_mram_array_u32(dpu, "curr_frontier", cf, lcf);
    dpu_insert_mram_array_u32(dpu, "node_levels", 0, lnl);

    // Copy COO data. Variable sized buffers must be copied last.
    uint32_t num_edges = coo[i].num_edges;
    dpu_set_u32(dpu, "num_edges", num_edges);
    dpu_insert_mram_array_u32(dpu, "nodes", coo[i].row_idxs, num_edges);
    dpu_insert_mram_array_u32(dpu, "neighbors", coo[i].col_idxs, num_edges);

    // Cache some MRAM addresses (address must be the same for all DPUs).
    DPU_ASSERT(dpu_copy_from(dpu, "next_frontier", 0, &nf_addr, sizeof(mram_addr_t)));
    DPU_ASSERT(dpu_copy_from(dpu, "curr_frontier", 0, &cf_addr, sizeof(mram_addr_t)));
  }

#if BENCHMARK_TIME
  t = clock() - t;
  pop_mram_time += ((double)t) / CLOCKS_PER_SEC;
#endif

  // Free resources.
  free(frontier);
  for (int i = 0; i < num_dpu; ++i)
    free_coo(coo[i]);

  // Start BFS algorithm.
  PRINT_INFO("Starting BFS algorithm.");
  if (prt == Row)
    bfs_vtx_row(len_cf, len_nf);
  else if (prt == Col)
    bfs_vtx_col(len_cf, len_nf);
  else
    bfs_vtx_2d(len_frontier, len_cf, len_nf, col_div);

  // Print node levels.
  print_node_levels(total_nodes, len_nl, 1);
}

// Cache DPU variable symbols for better performance.
void cache_symbols(struct dpu_program_t *program) {
  DPU_ASSERT(dpu_get_symbol(program, "__sys_used_mram_end", &mram_heap_sym));
  DPU_ASSERT(dpu_get_symbol(program, "level", &level_sym));
}

int main(int argc, char **argv) {

  enum Algorithm alg = SrcVtx;
  enum Partition prt = Row;
  char *bin_path;
  char *file = NULL;
  char *out_file = NULL;
  parse_args(argc, argv, &num_dpu, &alg, &prt, &bin_path, &file, &out_file);
  out = fopen(out_file, "w");

  PRINT_INFO("Allocating %u DPUs, %u tasklets each. Using %u bytes blocks for MRAM DMA.", num_dpu, NR_TASKLETS, BLOCK_SIZE);
  struct dpu_program_t *program;
  DPU_ASSERT(dpu_alloc(num_dpu, NULL, &set));
  DPU_ASSERT(dpu_load(set, bin_path, &program));
  cache_symbols(program);

  struct COO coo = load_coo(file, num_dpu);
  struct COO *coo_prts = partition_coo(coo, num_dpu, prt);
  free_coo(coo);

  if (alg == SrcVtx)
    start_src_vtx(coo_prts, num_dpu, prt);
  else if (alg == DstVtx) {
    start_dst_vtx(coo_prts, num_dpu, prt);
  } else if (alg == Edge) {
    start_edge(coo_prts, num_dpu, prt);
  }

  fclose(out);
  DPU_ASSERT(dpu_free(set));
  PRINT_INFO("Done");

#if BENCHMARK_TIME
  printf("dpu_compute_time %f\thost_comm_time %f\tpop_mram_time %f\tfetch_res_time %f\ttotal_alg %f\ttotal_pop_fetch %f\ttotal_all %f\n",
         dpu_compute_time, host_comm_time, pop_mram_time, fetch_res_time, dpu_compute_time + host_comm_time,
         pop_mram_time + fetch_res_time, dpu_compute_time + host_comm_time + pop_mram_time + fetch_res_time);
#endif

  return 0;
}
