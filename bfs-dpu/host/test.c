#include <dpu.h>
#include <dpu_log.h>
#include <dpu_memory.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#ifndef DPU_EXE
#define DPU_EXE "bin/testdpu"
#endif

void dpu_insert_mram_array_u32(struct dpu_set_t dpu, const char *symbol_name, uint32_t *src, uint32_t length) {
  // Get end of used MRAM pointer.
  mram_addr_t p_used_mram_end;
  DPU_ASSERT(dpu_copy_from(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Set the array pointer as the previous pointer.
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Copy the data to MRAM.
  size_t size = length * sizeof(uint32_t);
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_used_mram_end, (const uint8_t *)src, size));

  // Increment end of used MRAM pointer.
  p_used_mram_end += size;
  DPU_ASSERT(dpu_copy_to(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));
}

void dpu_insert_to_mram_uint64_array(struct dpu_set_t dpu, const char *symbol_name, uint64_t *src, uint64_t length) {
  // Get end of used MRAM pointer.
  mram_addr_t p_used_mram_end;
  DPU_ASSERT(dpu_copy_from(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Set the array pointer as the previous pointer.
  DPU_ASSERT(dpu_copy_to(dpu, symbol_name, 0, &p_used_mram_end, sizeof(mram_addr_t)));

  // Copy the data to MRAM.
  size_t size = length * sizeof(uint64_t);
  DPU_ASSERT(dpu_copy_to_mram(dpu.dpu, p_used_mram_end, (const uint8_t *)src, size));

  // Increment end of used MRAM pointer.
  p_used_mram_end += size;
  DPU_ASSERT(dpu_copy_to(dpu, "p_used_mram_end", 0, &p_used_mram_end, sizeof(mram_addr_t)));
}

int main() {
  struct dpu_set_t set, dpu;
  int dpuCount = 1;

  DPU_ASSERT(dpu_alloc(dpuCount, NULL, &set));
  DPU_ASSERT(dpu_load(set, DPU_EXE, NULL));

  uint64_t numCount = 2048;
  uint64_t *nums = malloc(numCount * sizeof(uint64_t));
  for (int i = 0; i < numCount; ++i) {
    nums[i] = (rand() % 100) + 1;
  }
  printf("\n\n");

  DPU_FOREACH(set, dpu) {
    // DPU_ASSERT(dpu_copy_to(dpu, "numCount", 0, &numCount64, sizeof(uint32_t)));
    dpu_insert_to_mram_uint64_array(dpu, "nums", nums, numCount);
  }

  DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
  DPU_FOREACH(set, dpu) {
    DPU_ASSERT(dpu_log_read(dpu, stdout));
  }

  printf("\n");
  DPU_ASSERT(dpu_free(set));
  return 0;
}