""" Determine, given the datafile and the expected bfs results,
    for each algorithm-partitioning pair the optimal number of tasklets
    per DPU and the optmial MRAM DMA block size.

    USAGE:  python3 bench_cyles.py <coo_data_file> <expected_results>

    We havedetermined consistently good performance with NR_TASKLETS=11
    and BLOCK_SIZE=256 across all algorithm combinations.
    Recommended as default.
"""

import subprocess
import filecmp
import math
import os
import sys
import logging
logging.basicConfig(filename='benchmark_cycles.log', level=logging.INFO)

datafile = sys.argv[1]
expected_node_levels = sys.argv[2]

num_dpus = 8

nr_tasklets = [x for x in range(1, 25)]
block_sizes = [2**x for x in range(3, 10)]  # [8, 16, 32, ... 512]

# (algorithm, partitioning) pairs
algs = [("src", "row"), ("src", "col"), ("src", "2d"),
        ("dst", "row"), ("dst", "col"), ("dst", "2d"),
        ("edge", "row"), ("edge", "col"), ("edge", "2d")]

configs = []


def bench_dpu_cycles(alg, prt, nr_tsk, block_size):
    """ Runs BFS on selected datafile with the specified configuration,
        and returns the total DPU cycles.
    """
    id_str = f"{alg}_{prt}_{nr_tsk}_{block_size}"
    res = f"res_{id_str}"

    make = f"NR_TASKLETS={nr_tsk} BLOCK_SIZE={block_size} BENCHMARK_CYCLES=true make all"
    run = f"./bin/bfs -n {num_dpus} -a {alg} -p {prt} -o {res} {datafile}"
    process = subprocess.run(make, shell=True)

    try:
        process = subprocess.run(
            run, shell=True, timeout=120, stdout=subprocess.PIPE, encoding="utf-8")
    except subprocess.TimeoutExpired:
        logging.error(f"BFS timeout ({id_str})")
        return

    if process.returncode > 0:
        logging.error(f"BFS failed to complete ({id_str})")
        return

    if not filecmp.cmp(res, expected_node_levels):
        logging.error(f"BFS output is incorrect ({id_str})")
        os.remove(res)
        return
    os.remove(res)

    total_cycles = 0
    for line in process.stdout.splitlines()[1:]:
        split_line = line.split(" ")
        total_cycles += int(split_line[0])

    return total_cycles


def get_best_config(alg, prt):
    """ Computes the optimal number of tasklets and block size
        for the specified (algorithm, partition) pair given the datafile.
    """

    least_cycles = math.inf
    best_nr_tasklets = 0
    best_block_size = 0

    for t in nr_tasklets:

        least_cycles_inner = math.inf
        best_block_size_inner = 0

        for b in block_sizes:
            dpu_cycles = bench_dpu_cycles(alg, prt, t, b)
            if dpu_cycles < least_cycles_inner:
                least_cycles_inner = dpu_cycles
                best_block_size_inner = b

        if dpu_cycles < least_cycles:
            least_cycles = dpu_cycles
            best_nr_tasklets = t
            best_block_size = best_block_size_inner

    return best_nr_tasklets, best_block_size


for a, p in algs:
    nr_tasklets, block_size = get_best_config(a, p)
    configs.append((a, p, nr_tasklets, block_size))

for a, p, t, b in configs:
    logging.info(
        f"Optimal config for {a}-{p} given datafile {datafile} is NR_TASKLETS={t} and BLOCK_SIZE={b}")
