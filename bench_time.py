
import subprocess
import logging
import filecmp
import os


def make(nr_tasklets, block_size):
    cmnd = f"NR_TASKLETS={nr_tasklets} BLOCK_SIZE={block_size} BENCHMARK_TIME=true make all"
    subprocess.run(cmnd, stdout=subprocess.PIPE, shell=True)


def get_max_degree(datafile):
    """
        Returns the node with the highest number of edges and its number of edges.
    """
    with open(datafile) as f:
        next(f)
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

    return max_degree_node, max_degree


def bfs(datafile, expected_node_levels, alg, prt, num_dpus):

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


success, dpu_compute_time, host_comm_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all = bfs(
    "data/loc-brightkite_edges.txt", "results/true/loc-brightkite_edges.txt", "src", "row", 64)
if success:
    print(f"dpu_compute_time={dpu_compute_time}, host_comm_time={host_comm_time}, pop_mram_time={pop_mram_time}, fetch_res_time={fetch_res_time}, total_alg={total_alg}, total_pop_fetch={total_pop_fetch}, total_all={total_all}")
