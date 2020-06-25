#include <alloc.h>
#include <barrier.h>
#include <defs.h>
#include <mram.h>
#include <mutex.h>
#include <perfcounter.h>
#include <seqread.h>
#include <stdint.h>
#include <stdio.h>

#define TILE_SIZE 512

__host __mram_ptr void *p_used_mram_end = DPU_MRAM_HEAP_POINTER; // Points to the end of used MRAM addresses.

__host __mram_ptr uint32_t *nums;

__host __dma_aligned uint32_t factorial = 1;

int main() {

  (void)perfcounter_config(COUNT_CYCLES, true);

  for (uint32_t i = 10; i > 0; i--)
    factorial *= i;

  mram_write(&factorial, (__mram_ptr void *)0, 8);

  perfcounter_t run_time = perfcounter_get();
  printf("\nRes: %u -- Perf: %lu ", factorial, run_time);

  return 0;
}