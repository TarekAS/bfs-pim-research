
#include <alloc.h>
#include <defs.h>
#include <ktrace.h>
#include <mram.h>
#include <seqread.h>

#define printf ktrace
#define PRINT_ERROR(fmt, ...) printf("\033[0;31mERROR:\033[0m   "fmt"\n", ##__VA_ARGS__)

#define ROUND_UP_TO_MULTIPLE_OF_2(x)    ((((x) - 1)/2 + 1)*2)
#define ROUND_UP_TO_MULTIPLE_OF_8(x)    ((((x) - 1)/8 + 1)*8)
#define MIN(x, y)   (((x) < (y))?(x):(y))

#define NUM_TASKLETS 16
#define TASKLETS_INITIALIZER \
    TASKLETS(NUM_TASKLETS, main, 1024, 0)
#include <rt.h>

struct Nonzero {
    uint32_t col;
    float value;
};

int main() {

    // Load parameters
    uint32_t numParams = 4;
    uint32_t paramsSize = ROUND_UP_TO_MULTIPLE_OF_8(numParams*sizeof(uint32_t)); /* 16 */
    mram_addr_t params_m = 0;
    uint32_t* params_w = (uint32_t*) mem_alloc_dma(paramsSize);
    mram_read16(params_m, params_w);
    uint32_t numRows = params_w[0]; /* Number of rows assigned to this DPU */
    uint32_t numCols = params_w[1]; /* Number of columns in the matrix */
    uint32_t rowPtrsOffset = params_w[2];  /* Offset of the row pointers */
    uint32_t numNonzeros = params_w[3];  /* Total number of nonzeros in the rows assigned to this DPU */

    // Sanity check
    if(me() == 0) {
        if(numRows%2 != 0) {
            // The number of rows assigned to the DPU must be a multiple of two to ensure that writes to the output vector are aligned to 8 bytes
            PRINT_ERROR("The number of rows is not a multiple of two!");
        }
    }

    // Identify tasklet's rows
    uint32_t numRowsPerTasklet = ROUND_UP_TO_MULTIPLE_OF_2((numRows - 1)/NUM_TASKLETS + 1); // Multiple of two to ensure that access to rowPtrs and outVector is 8-byte aligned
    uint32_t taskletRowsStart = me()*numRowsPerTasklet;
    uint32_t taskletNumRows;
    if(taskletRowsStart > numRows) {
        taskletNumRows = 0;
    } else if(taskletRowsStart + numRowsPerTasklet > numRows) {
        taskletNumRows = numRows - taskletRowsStart;
    } else {
        taskletNumRows = numRowsPerTasklet;
    }

    // Only process tasklets with nonzero number of rows
    if(taskletNumRows > 0) {

        // Extract CSR array pointers (all 8-byte aligned)
        mram_addr_t rowPtrs_m = params_m + paramsSize;
        uint32_t rowPtrsSize = ROUND_UP_TO_MULTIPLE_OF_8((numRows + 1)*sizeof(uint32_t));
        mram_addr_t nonzeros_m = rowPtrs_m + rowPtrsSize;
        uint32_t nonzerosSize = ROUND_UP_TO_MULTIPLE_OF_8(numNonzeros*sizeof(struct Nonzero));
        mram_addr_t inVector_m = nonzeros_m + nonzerosSize;
        uint32_t inVectorSize = ROUND_UP_TO_MULTIPLE_OF_8(numCols*sizeof(float));
        mram_addr_t outVector_m = inVector_m + inVectorSize;
        uint32_t outVectorSize = ROUND_UP_TO_MULTIPLE_OF_8(numRows*sizeof(float));

        // Initialize row pointer sequential reader
        mram_addr_t taskletRowPtrs_m = rowPtrs_m + taskletRowsStart*sizeof(uint32_t);
        uint32_t* taskletRowPtrs_w;
        seqreader_t rowPtrReader = seqread_init(seqread_alloc(), taskletRowPtrs_m, (void**) &taskletRowPtrs_w);
        uint32_t firstRowPtr = *taskletRowPtrs_w;

        // Initialize nonzeros sequential reader
        uint32_t taskletNonzerosStart = firstRowPtr - rowPtrsOffset;
        mram_addr_t taskletNonzeros_m = nonzeros_m + taskletNonzerosStart*sizeof(struct Nonzero); // 8-byte aligned because Nonzero is 8 bytes
        struct Nonzero* taskletNonzeros_w;
        seqreader_t nonzerosReader = seqread_init(seqread_alloc(), taskletNonzeros_m, (void**) &taskletNonzeros_w);

        // Initialize input vector cache
        uint32_t inVectorTileSize = 64;
        float* inVectorTile_w = mem_alloc_dma(inVectorTileSize*sizeof(float));
        mram_read256(inVector_m, inVectorTile_w);
        uint32_t currInVectorTileIdx = 0;

        // Initialize output vector cache
        mram_addr_t taskletOutVector_m = outVector_m + taskletRowsStart*sizeof(float);
        uint32_t outVectorTileSize = 64;
        float* outVectorTile_w = mem_alloc_dma(outVectorTileSize*sizeof(float));

        // SpMV
        uint32_t nextRowPtr = firstRowPtr;
        for(uint32_t row = 0; row < taskletNumRows; ++row) {

            // Find row nonzeros
            taskletRowPtrs_w = seqread_get(taskletRowPtrs_w, sizeof(uint32_t), &rowPtrReader);
            uint32_t rowPtr = nextRowPtr;
            nextRowPtr = *taskletRowPtrs_w;
            uint32_t taskletNNZ = nextRowPtr - rowPtr;

            // Multiply row with vector
            float outValue = 0.0f;
            for(uint32_t nzIdx = 0; nzIdx < taskletNNZ; ++nzIdx) {

                // Get matrix value
                float matValue = taskletNonzeros_w->value;

                // Get input vector value
                uint32_t col = taskletNonzeros_w->col;
                uint32_t inVectorTileIdx = col/inVectorTileSize;
                uint32_t inVectorTileOffset = col%inVectorTileSize;
                if(inVectorTileIdx != currInVectorTileIdx) {
                    mram_read256(inVector_m + inVectorTileIdx*inVectorTileSize*sizeof(float), inVectorTile_w);
                    currInVectorTileIdx = inVectorTileIdx;
                }
                float inValue = inVectorTile_w[inVectorTileOffset];

                // Multiply and add
                outValue += matValue*inValue;

                // Read next nonzero
                taskletNonzeros_w = seqread_get(taskletNonzeros_w, sizeof(struct Nonzero), &nonzerosReader); // Last read will be out of bounds and unused

            }

            // Store output
            uint32_t outVectorTileIdx = row/outVectorTileSize;
            uint32_t outVectorTileOffset = row%outVectorTileSize;
            outVectorTile_w[outVectorTileOffset] = outValue;
            if(outVectorTileOffset == outVectorTileSize - 1) { // Last element in tile
                mram_write256(outVectorTile_w, taskletOutVector_m + outVectorTileIdx*outVectorTileSize*sizeof(float));
            } else if(row == taskletNumRows - 1) { // Last row for tasklet
                mram_writeX(outVectorTile_w, taskletOutVector_m + outVectorTileIdx*outVectorTileSize*sizeof(float), (taskletNumRows%outVectorTileSize)*sizeof(float));
            }

        }
    }

    return 0;

}

