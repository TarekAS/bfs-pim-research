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


def bench_dpu_cycles(alg, prt, nr_tsk, block_size):
    """ Runs BFS on selected datafile with the specified configuration,
        and returns the total DPU cycles.
    """
    id_str = f"{alg}_{prt}_{nr_tsk}_{block_size}"
    res = f"res_{id_str}"

    try:
        run = f"./bin/bfs -n {num_dpus} -a {alg} -p {prt} -o {res} {datafile}"
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


datafile = sys.argv[1]
expected_node_levels = sys.argv[2]
logging.basicConfig(
    filename='bench_cycles.log',
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

num_dpus = 64
nr_tasklets = [x for x in range(1, 25)]
block_sizes = [2**x for x in range(3, 10)]  # [8, 16, 32, ... 512]

# (algorithm, partitioning) pairs
algs = [("src", "row"), ("src", "col"), ("src", "2d"),
        ("dst", "row"), ("dst", "col"), ("dst", "2d"),
        ("edge", "row"), ("edge", "col"), ("edge", "2d")]

configs = {}


for t in nr_tasklets:

    for b in block_sizes:
        make = f"NR_TASKLETS={t} BLOCK_SIZE={b} BENCHMARK_CYCLES=true make all"
        process = subprocess.run(make, shell=True)

        for a, p in algs:
            dpu_cycles = bench_dpu_cycles(a, p, t, b)
            print(f"{a} {p} t={t} b={b} {dpu_cycles}")
