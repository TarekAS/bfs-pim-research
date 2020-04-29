
#include <dpu.h>
#include <dpu_log.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define DPU_BINARY "./bin/task.bin"

#define DPU_LOG_DIRECTORY "./log/host/"

#define PRINT_ERROR(fmt, ...)   fprintf(stderr, "\033[0;31mERROR:\033[0m   "fmt"\n", ##__VA_ARGS__)
#define PRINT_WARNING(fmt, ...) fprintf(stderr, "\033[0;35mWARNING:\033[0m "fmt"\n", ##__VA_ARGS__)
#define PRINT_INFO(fmt, ...)    fprintf(stderr, "\033[0;32mINFO:\033[0m    "fmt"\n", ##__VA_ARGS__)

#ifdef DPU_PROFILE
#define STR(a) #a
#define STRINGIFY(a) STR(a)
#define DPU_PROFILE_STR STRINGIFY(DPU_PROFILE)
#else
#define DPU_PROFILE_STR ""
#endif

#define ROUND_UP_TO_MULTIPLE_OF_2(x)    ((((x) + 1)/2)*2)
#define ROUND_UP_TO_MULTIPLE_OF_8(x)    ((((x) + 7)/8)*8)

#define DPU_CAPACITY (64 << 20) // A DPU's capacity is 64 MiB

struct Nonzero {
    uint32_t col;
    float value;
};

struct COOMatrix {
    uint32_t numRows;
    uint32_t numCols;
    uint32_t numNonzeros;
    uint32_t* rowIdxs;
    struct Nonzero* nonzeros;
};

struct CSRMatrix {
    uint32_t numRows;
    uint32_t numCols;
    uint32_t numNonzeros;
    uint32_t* rowPtrs;
    struct Nonzero* nonzeros;
};

struct COOMatrix readCOOMatrix(const char* fileName) {

    struct COOMatrix cooMatrix;

    // Initialize fields
    FILE* fp = fopen(fileName, "r");
    fscanf(fp, "%u", &cooMatrix.numRows);
    if(cooMatrix.numRows%2 == 1) {
        PRINT_WARNING("Reading matrix %s: number of rows must be even. Padding with an extra row.", fileName);
        cooMatrix.numRows++;
    }
    fscanf(fp, "%u", &cooMatrix.numCols);
    fscanf(fp, "%u", &cooMatrix.numNonzeros);
    cooMatrix.rowIdxs = malloc(cooMatrix.numNonzeros*sizeof(uint32_t));
    cooMatrix.nonzeros = malloc(cooMatrix.numNonzeros*sizeof(struct Nonzero));

    PRINT_INFO("Reading matrix %s: %u rows, %u columns, %u nonzeros", fileName, cooMatrix.numRows, cooMatrix.numCols, cooMatrix.numNonzeros);

    // Read the nonzeros
    for(uint32_t i = 0; i < cooMatrix.numNonzeros; ++i) {
        uint32_t rowIdx;
        fscanf(fp, "%u", &rowIdx);
        cooMatrix.rowIdxs[i] = rowIdx - 1; // File format indexes begin at 1
        uint32_t colIdx;
        fscanf(fp, "%u", &colIdx);
        cooMatrix.nonzeros[i].col = colIdx - 1; // File format indexes begin at 1
        fscanf(fp, "%f", &cooMatrix.nonzeros[i].value);
    }

    return cooMatrix;

}

void freeCOOMatrix(struct COOMatrix cooMatrix) {
    free(cooMatrix.rowIdxs);
    free(cooMatrix.nonzeros);
}

struct CSRMatrix coo2csr(struct COOMatrix cooMatrix) {

    struct CSRMatrix csrMatrix;

    // Initialize fields
    csrMatrix.numRows = cooMatrix.numRows;
    csrMatrix.numCols = cooMatrix.numCols;
    csrMatrix.numNonzeros = cooMatrix.numNonzeros;
    csrMatrix.rowPtrs = malloc((csrMatrix.numRows + 1)*sizeof(uint32_t));
    csrMatrix.nonzeros = malloc(csrMatrix.numNonzeros*sizeof(struct Nonzero));

    // Histogram rowIdxs
    memset(csrMatrix.rowPtrs, 0, (csrMatrix.numRows + 1)*sizeof(uint32_t));
    for(uint32_t i = 0; i < cooMatrix.numNonzeros; ++i) {
        uint32_t rowIdx = cooMatrix.rowIdxs[i];
        csrMatrix.rowPtrs[rowIdx]++;
    }

    // Prefix sum rowPtrs
    uint32_t sumBeforeNextRow = 0;
    for(uint32_t rowIdx = 0; rowIdx < csrMatrix.numRows; ++rowIdx) {
        uint32_t sumBeforeRow = sumBeforeNextRow;
        sumBeforeNextRow += csrMatrix.rowPtrs[rowIdx];
        csrMatrix.rowPtrs[rowIdx] = sumBeforeRow;
    }
    csrMatrix.rowPtrs[csrMatrix.numRows] = sumBeforeNextRow;

    // Bin the nonzeros
    for(uint32_t i = 0; i < cooMatrix.numNonzeros; ++i) {
        uint32_t rowIdx = cooMatrix.rowIdxs[i];
        uint32_t nnzIdx = csrMatrix.rowPtrs[rowIdx]++;
        csrMatrix.nonzeros[nnzIdx] = cooMatrix.nonzeros[i];
    }

    // Restore rowPtrs
    for(uint32_t rowIdx = csrMatrix.numRows - 1; rowIdx > 0; --rowIdx) {
        csrMatrix.rowPtrs[rowIdx] = csrMatrix.rowPtrs[rowIdx - 1];
    }
    csrMatrix.rowPtrs[0] = 0;

    return csrMatrix;

}

void freeCSRMatrix(struct CSRMatrix csrMatrix) {
    free(csrMatrix.rowPtrs);
    free(csrMatrix.nonzeros);
}

void initVector(float* vec, uint32_t size) {
    for(uint32_t i = 0; i < size; ++i) {
        vec[i] = 1.0f;
    }
}

void usage() {
    PRINT_ERROR(
            "\nUsage:  ./program [options]"
            "\n"
            "\nGeneral options:"
            "\n    -h        help"
            "\n    -d <D>    DPU type (default=fsim)"
            "\n    -r <R>    # of ranks (default=8)"
            "\n"
            "\nBenchmark-specific options:"
            "\n    -f <F>    input matrix file name (default=data/bcsstk13.mtx)"
            "\n");
}

typedef struct Params {
  char* dpuType;
  int   numRanks;
  int   n_warmup;
  int   n_reps;
  char* fileName;
} Params;

struct Params input_params(int argc, char **argv) {
    struct Params p;
    p.dpuType       = "fsim";
    p.numRanks      = 8;
    p.n_warmup      = 2;
    p.n_reps        = 5;
    p.fileName      = "data/bcsstk13.mtx";

    int opt;
    while((opt = getopt(argc, argv, "hd:r:i:")) >= 0) {
        switch(opt) {
            case 'h':
                      usage();
                      exit(0);
            case 'd': p.dpuType       = optarg; break;
            case 'r': p.numRanks      = atoi(optarg); break;
            case 'f': p.fileName      = optarg; break;
            default:
                      PRINT_ERROR("Unrecognized option!");
                      usage();
                      exit(0);
        }
    }
    assert(p.numRanks > 0 && "Invalid # of ranks!");

    return p;
}

void allocateDPUs(struct dpu_rank_t **rank, dpu_type_t dpuType) {

    struct dpu_logging_config_t log_config = {
        .source = KTRACE,
        .destination_directory_name = DPU_LOG_DIRECTORY
    };

    struct dpu_param_t params = {
        .type = dpuType,
        .profile = DPU_PROFILE_STR,
        .logging_config = &log_config
    };

    if (dpu_alloc(&params, rank) != DPU_API_SUCCESS) {
        PRINT_ERROR("Cannot allocate rank");
        exit(0);
    }

    if (dpu_load_all(*rank, DPU_BINARY) != DPU_API_SUCCESS) {
        dpu_free(*rank);
        PRINT_ERROR("Cannot load DPU program");
        exit(0);
    }

}

void copyToDPU(struct dpu_t* dpu, uint8_t* hostPtr, uint32_t mramIdx, uint32_t size) {
    if(dpu_copy_to_dpu(dpu, hostPtr, mramIdx, size) != DPU_API_SUCCESS) {
        PRINT_ERROR("Cannot copy to DPU at MRAM index %u, size %u", mramIdx, size);
        exit(0);
    }
}

void copyFromDPU(struct dpu_t* dpu, uint32_t mramIdx, uint8_t* hostPtr, uint32_t size) {
    if(dpu_copy_from_dpu(dpu, mramIdx, hostPtr, size) != DPU_API_SUCCESS) {
        PRINT_ERROR("Cannot copy from DPU");
        exit(0);
    }
}

int main(int argc, char** argv) {

    // Process parameters
    struct Params p = input_params(argc, argv);

    // DPU type
    dpu_type_t dpuType;
    if (strcmp(p.dpuType, "fsim") == 0) {
        dpuType = FUNCTIONAL_SIMULATOR;
        PRINT_INFO("DPU type: functional simulator");
    } else if (strcmp(p.dpuType, "asic") == 0) {
        dpuType = HW;
        PRINT_INFO("DPU type: ASIC");
    } else if (strcmp(p.dpuType, "fpga") == 0) {
        dpuType = HW;
        PRINT_INFO("DPU type: FPGA");
    }

    // Allocate DPUs
    uint32_t numRanks = p.numRanks;
    uint32_t numDPUsPerRank;
    struct dpu_rank_t* dpuRanks[numRanks];
    PRINT_INFO("Allocating %u ranks", numRanks);
    for(uint32_t r = 0; r < numRanks; ++r) {
        allocateDPUs(&dpuRanks[r], dpuType);
        if(dpu_get_nr_of_dpus_in(dpuRanks[0], &numDPUsPerRank) != DPU_API_SUCCESS) {
            PRINT_ERROR("Could not get number of DPUs in rank");
            exit(0);
        } else {
            PRINT_INFO("    Allocated %u DPU(s) in rank %u.", numDPUsPerRank, r);
        }
    }

    // Initialize SpMV data structures
    struct COOMatrix cooMatrix = readCOOMatrix(p.fileName);
    struct CSRMatrix csrMatrix = coo2csr(cooMatrix);
    uint32_t numRows = csrMatrix.numRows;
    uint32_t numCols = csrMatrix.numCols;
    uint32_t numNonzeros = csrMatrix.numNonzeros;
    uint32_t* rowPtrs = csrMatrix.rowPtrs;
    struct Nonzero* nonzeros = csrMatrix.nonzeros;
    float* inVector = malloc(csrMatrix.numCols*sizeof(float));
    initVector(inVector, csrMatrix.numCols);
    float* outVector = malloc(csrMatrix.numRows*sizeof(float));

    // Partition data structure across DPUs
    uint32_t numDPUs = numRanks*numDPUsPerRank;
    uint32_t numRowsPerDPU = ROUND_UP_TO_MULTIPLE_OF_2((numRows - 1)/numDPUs + 1);
    PRINT_INFO("Assigning %u rows per DPU", numRowsPerDPU);
    uint32_t dpuOutVector_m[numDPUs];
    for(uint32_t r = 0; r < numRanks; ++r) {
        struct dpu_t *dpu;
        uint32_t d = 0;
        DPU_FOREACH(dpuRanks[r], dpu) {
            uint32_t dpuIdx = r*numDPUsPerRank + d;
            uint32_t dpuStartRowIdx = dpuIdx*numRowsPerDPU;
            uint32_t dpuNumRows;
            if(dpuStartRowIdx > numRows) {
                dpuNumRows = 0;
            } else if(dpuStartRowIdx + numRowsPerDPU > numRows) {
                dpuNumRows = numRows - dpuStartRowIdx;
            } else {
                dpuNumRows = numRowsPerDPU;
            }
            PRINT_INFO("    DPU %u (rank %u, DPU %u):", dpuIdx, r, d);
            PRINT_INFO("        Receives %u rows", dpuNumRows);
            if(dpuNumRows > 0) {
                // Find DPU's CSR matrix partition
                uint32_t* dpuRowPtrs_h = rowPtrs + dpuStartRowIdx;
                uint32_t dpuRowPtrsOffset = dpuRowPtrs_h[0];
                struct Nonzero* dpuNonzeros_h = nonzeros + dpuRowPtrsOffset;
                uint32_t dpuNumNonzeros = dpuRowPtrs_h[dpuNumRows] - dpuRowPtrsOffset;
                // Layout DPU memory
                uint32_t dpuParams_m = 0;
                uint32_t numParams = 4;
                uint32_t paramsSize = ROUND_UP_TO_MULTIPLE_OF_8(numParams*sizeof(uint32_t)); /* 16 */
                uint32_t dpuRowPtrs_m = dpuParams_m + paramsSize;
                uint32_t dpuRowPtrsSize = (dpuNumRows + 1)*sizeof(uint32_t);
                uint32_t dpuRowPtrsSizePadded = ROUND_UP_TO_MULTIPLE_OF_8(dpuRowPtrsSize);
                uint32_t dpuNonzeros_m = dpuRowPtrs_m + dpuRowPtrsSizePadded;
                uint32_t dpuNonzerosSize = dpuNumNonzeros*sizeof(struct Nonzero);
                uint32_t dpuNonzerosSizePadded = ROUND_UP_TO_MULTIPLE_OF_8(dpuNonzerosSize);
                uint32_t dpuInVector_m = dpuNonzeros_m + dpuNonzerosSizePadded;
                uint32_t dpuInVectorSize = numCols*sizeof(float);
                uint32_t dpuInVectorSizePadded = ROUND_UP_TO_MULTIPLE_OF_8(dpuInVectorSize);
                dpuOutVector_m[dpuIdx] = dpuInVector_m + dpuInVectorSizePadded;
                uint32_t dpuOutVectorSize = dpuNumRows*sizeof(float);
                if(dpuOutVectorSize%8 != 0) {
                    PRINT_WARNING("        Output sub-vector is not a multiple of 8 bytes!");
                }
                uint32_t totalDPUMemory = dpuOutVector_m[dpuIdx] + dpuOutVectorSize;
                if(totalDPUMemory <= DPU_CAPACITY) {
                    PRINT_INFO("        Total memory allocated is %d bytes", totalDPUMemory);
                } else {
                    PRINT_ERROR("        Total memory allocated is %d bytes which exceeds the DPU capacity (%d bytes)!", totalDPUMemory, DPU_CAPACITY);
                    exit(0);
                }
                // Send data to DPU
                PRINT_INFO("        Copying data to DPU");
                copyToDPU(dpu, (uint8_t*)&dpuNumRows, dpuParams_m + 0, sizeof(uint32_t));
                copyToDPU(dpu, (uint8_t*)&numCols, dpuParams_m + sizeof(uint32_t), sizeof(uint32_t));
                copyToDPU(dpu, (uint8_t*)&dpuRowPtrsOffset, dpuParams_m + 2*sizeof(uint32_t), sizeof(uint32_t));
                copyToDPU(dpu, (uint8_t*)&dpuNumNonzeros, dpuParams_m + 3*sizeof(uint32_t), sizeof(uint32_t));
                copyToDPU(dpu, (uint8_t*)dpuRowPtrs_h, dpuRowPtrs_m, dpuRowPtrsSize);
                copyToDPU(dpu, (uint8_t*)dpuNonzeros_h, dpuNonzeros_m, dpuNonzerosSize);
                copyToDPU(dpu, (uint8_t*)inVector, dpuInVector_m, dpuInVectorSize);
            }
            ++d;
        }
    }

    // Boot all DPUs in all ranks
    PRINT_INFO("Booting DPUs");
    for(uint32_t r = 0; r < numRanks; ++r) {
        dpu_boot_all(dpuRanks[r], ASYNCHRONOUS);
    }

    // Wait for all DPUs to finish
    PRINT_INFO("Waiting for DPUs to finish");
    for(uint32_t r = 0; r < numRanks; ++r) {
        uint32_t numDPUsRunning = 0;
        dpu_run_status_t *status = calloc(sizeof(dpu_run_status_t), numDPUsPerRank);
        do {
            dpu_get_all_status(dpuRanks[r], status, &numDPUsRunning);
        } while (numDPUsRunning != 0);
        free(status);
    }

    // Copy back result
    PRINT_INFO("Copying back the result");
    for(uint32_t r = 0; r < numRanks; ++r) {
        struct dpu_t *dpu;
        uint32_t d = 0;
        DPU_FOREACH(dpuRanks[r], dpu) {
            uint32_t dpuIdx = r*numDPUsPerRank + d;
            uint32_t dpuStartRowIdx = dpuIdx*numRowsPerDPU;
            uint32_t dpuNumRows;
            if(dpuStartRowIdx > numRows) {
                dpuNumRows = 0;
            } else if(dpuStartRowIdx + numRowsPerDPU > numRows) {
                dpuNumRows = numRows - dpuStartRowIdx;
            } else {
                dpuNumRows = numRowsPerDPU;
            }
            if(dpuNumRows > 0) {
                copyFromDPU(dpu, dpuOutVector_m[dpuIdx], (uint8_t*)(outVector + dpuStartRowIdx), dpuNumRows*sizeof(float));
            }
            ++d;
        }
    }

    // Verify the result
    PRINT_INFO("Verifying the result");
    for(uint32_t rowIdx = 0; rowIdx < numRows; ++rowIdx) {
        float sum = 0.0f;
        for(uint32_t i = rowPtrs[rowIdx]; i < rowPtrs[rowIdx + 1]; ++i) {
            uint32_t colIdx = nonzeros[i].col;
            float value = nonzeros[i].value;
            sum += inVector[colIdx]*value;
        }
        float diff = (sum - outVector[rowIdx])/outVector[rowIdx];
        const float tolerance = 0.00001;
        if(diff > tolerance || diff < -tolerance) {
            PRINT_ERROR("Mismatch at index %u (CPU result = %f, DPU result = %f)", rowIdx, sum, outVector[rowIdx]);
        }
    }

    // Display DPU Logs
    printf("Displaying DPU Logs:\n");
    for (unsigned int r = 0; r < numRanks; ++r) {
        struct dpu_t *dpu;
        unsigned int d = 0;
        DPU_FOREACH(dpuRanks[r], dpu) {
            printf("Rank %d, DPU %d:\n", r, d);
            if(!dpulog_read_for_dpu(dpu, stdout)) {
                PRINT_ERROR("cannot display DPU log correctly");
            }
        }
        ++d;
    }

    // Deallocate data structures
    freeCOOMatrix(cooMatrix);
    freeCSRMatrix(csrMatrix);
    free(inVector);
    free(outVector);

    return 0;

}

