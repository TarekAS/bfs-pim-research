"""
    Benchmarks bfs-dpu.
    Usage:
        For a single datafile:  $ python3 bench_time.py <datafile> <expected_node_levels>
        For multiple datafiles: $ python3 bench_time.py
            Put the datafiles in ./data and the expected_node_levels in ./results/cpu.
    
    The expected_node_levels are used to verify the correctness of the BFS output. 

    Prints timing results to bench_time_results.
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
        return False, 0, 0, 0, 0, 0, 0, 0

    if process.returncode > 0:
        logging.error(f"BFS failed to complete ({id_str})")
        if os.path.exists(res):
            os.remove(res)
        return False, 0, 0, 0, 0, 0, 0, 0

    if not filecmp.cmp(res, expected_node_levels):
        logging.error(f"BFS output is incorrect ({id_str})")
        if os.path.exists(res):
            os.remove(res)
        return False, 0, 0, 0, 0, 0, 0, 0

    if os.path.exists(res):
        os.remove(res)

    times = process.stdout.split("\t")

    dpu_compute_time = float(times[0].split(" ")[1])
    host_comm_time = float(times[1].split(" ")[1])
    pop_mram_time = float(times[2].split(" ")[1])
    fetch_res_time = float(times[3].split(" ")[1])
    total_alg = float(times[4].split(" ")[1])
    total_pop_fetch = float(times[5].split(" ")[1])
    total_all = float(times[6].split(" ")[1])

    return True, dpu_compute_time, host_comm_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all


logging.basicConfig(filename='bench_time.error.log', level=logging.ERROR)
make(11, 32)


# (algorithm, partitioning) pairs
algs = [("src", "row"), ("src", "col"), ("src", "2d"),
        ("dst", "row"), ("dst", "col"), ("dst", "2d"),
        ("edge", "row"), ("edge", "col"), ("edge", "2d")]

datafiles = []
if len(sys.argv) > 1:
    datafiles = [(sys.argv[1], sys.argv[2])]
else:
    data = os.listdir("data")
    expected = os.listdir("results/cpu")
    missing = set(data) - set(expected)
    if len(missing) > 0:
        for f in missing:
            print(f"Missing expected output of file {f}")
        exit()
    for f in data:
        datafiles.append((f"data/{f}", f"results/cpu/{f}"))

min_dpus = 8
max_dpus = 640

outfile = "bench_time_results"
counter = 2
while os.path.isfile(outfile):
    outfile = f"bench_time_results_{counter}"
    counter += 1

f = open(outfile, "w+")

f.write("datafile\tsuccess\tnum_nodes\tnum_edges\tmax_degree_node\tmax_degree\tdpu_compute_time\thost_comm_time\tpop_mram_time\tfetch_res_time\ttotal_alg\ttotal_pop_fetch\ttotal_all")
f.flush()

for datafile, expected in datafiles:
    for alg, prt in algs:
        for num_dpus in range(min_dpus, max_dpus+8, 8):

            success, dpu_compute_time, host_comm_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all = bfs(
                datafile, expected, alg, prt, num_dpus)
            num_nodes, num_edges, max_degree_node, max_degree = get_metadata(
                datafile)

            f.write(f"{os.path.basename(datafile)}\t{alg}\t{prt}\t{num_dpus}\t{success}\t{num_nodes}\t{num_edges}\t{max_degree_node}\t{max_degree}\t{dpu_compute_time}\t{host_comm_time}\t{pop_mram_time}\t{fetch_res_time}\t{total_alg}\t{total_pop_fetch}\t{total_all}")
            f.flush()
f.close()
