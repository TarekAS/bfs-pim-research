"""
    Benchmarks bfs-dpu.
    Usage:
        For a single datafile:  $ python3 bench_time.py <datafile> <expected_node_levels>
        For multiple datafiles: $ python3 bench_time.py
            Put the datafiles in ./data and the expected_node_levels in ./results/expected.

    The expected_node_levels are used to verify the correctness of the BFS output.

    Prints timing results to bench_results.
"""


import subprocess
import logging
import filecmp
import sys
import os


def make(nr_tasklets, block_size):
    cmnd = f"NR_TASKLETS={nr_tasklets} BLOCK_SIZE={block_size} BENCHMARK_TIME=true make all"
    subprocess.run(cmnd, stdout=subprocess.PIPE, shell=True)


def get_metadata(datafile):
    with open(datafile) as f:
        num_nodes, num_edges = next(f).split()
        max_degree = 0
        max_degree_node = 0
        current_degree = 0
        current_node = 0

        for line in f:
            node = line.split()[0]

            if (node != current_node):
                current_node = node
                current_degree = 0

            current_degree += 1
            if current_degree > max_degree:
                max_degree = current_degree
                max_degree_node = node

    return num_nodes, num_edges, max_degree_node, max_degree


def bfs(datafile, expected_node_levels, alg, prt, num_dpus):
    """
        Runs specified BFS algorithm on a datafile using the specified number
        of DPUs. The output is compared with the expected_node_levels to verify
        that the run was correct.
    """

    id_str = f"{alg}_{prt}_{os.path.basename(datafile)}_{num_dpus}"
    res = f"res_{id_str}"

    run = f"./bin/bfs -n {num_dpus} -a {alg} -p {prt} -o {res} {datafile}"

    try:
        process = subprocess.run(
            run, shell=True, stdout=subprocess.PIPE, encoding="utf-8")
    except Exception:
        logging.error(f"BFS failed to run ({id_str})")
        return False, 0, 0, 0, 0, 0, 0, 0, 0

    if process.returncode > 0:
        logging.error(f"BFS failed to complete ({id_str})")
        if os.path.exists(res):
            os.remove(res)
        return False, 0, 0, 0, 0, 0, 0, 0, 0

    times = process.stdout.split(" ")

    dpu_compute_time = float(times[1])
    host_comm_time = float(times[3])
    host_aggr_time = float(times[5])
    pop_mram_time = float(times[7])
    fetch_res_time = float(times[9])
    total_alg = float(times[11])
    total_pop_fetch = float(times[13])
    total_all = float(times[15])

    if not filecmp.cmp(res, expected_node_levels):
        logging.error(f"BFS output is incorrect ({id_str})")
        if os.path.exists(res):
            os.remove(res)
        return False, dpu_compute_time, host_comm_time, host_aggr_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all

    if os.path.exists(res):
        os.remove(res)

    return True, dpu_compute_time, host_comm_time, host_aggr_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all


logging.basicConfig(filename='bench.error.log', level=logging.ERROR)

# Compile code with 11 Tasklets and 32 bytes block size.
make(11, 32)

# (algorithm, partitioning) pairs
algs = [("top", "row"), ("top", "col"), ("top", "2d"),
        ("bot", "row"), ("bot", "col"), ("bot", "2d"),
        ("edge", "row"), ("edge", "col"), ("edge", "2d")]

# Get datafiles from args.
datafiles = []
if len(sys.argv) > 1:
    datafiles = [(sys.argv[1], sys.argv[2])]
else:
    data = os.listdir("data")
    expected = os.listdir("results/expected")
    missing = set(data) - set(expected)
    if len(missing) > 0:
        for f in missing:
            print(f"Missing expected output of file {f}")
        exit()
    for f in data:
        datafiles.append((f"data/{f}", f"results/expected/{f}"))

# Create unique output file.
outfile = "bench_results"
counter = 2
while os.path.isfile(outfile):
    outfile = f"bench_results_{counter}"
    counter += 1
f = open(outfile, "w+")

# Write header.
f.write("success datafile alg prt num_dpus dpu_compute_time host_comm_time host_aggr_time pop_mram_time fetch_res_time total_alg total_pop_fetch total_all\n")
f.flush()

# Run benchmarks on each datafile, for each combination of bfs variation and dpu count.
dpu_count = [8, 16, 32, 64, 128, 256]
for datafile, expected in datafiles:
    for alg, prt in algs:
        for num_dpus in dpu_count:

            success, dpu_compute_time, host_comm_time, host_aggr_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all = bfs(
                datafile, expected, alg, prt, num_dpus)
            # num_nodes, num_edges, max_degree_node, max_degree = get_metadata(datafile)

            f.write(f"{success} {os.path.basename(datafile)} {alg} {prt} {num_dpus} {dpu_compute_time} {host_comm_time} {host_aggr_time} {pop_mram_time} {fetch_res_time} {total_alg} {total_pop_fetch} {total_all}\n")
            f.flush()
f.close()
