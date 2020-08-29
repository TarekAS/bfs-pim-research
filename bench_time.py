
import subprocess
import logging
import filecmp
import os
logging.basicConfig(filename='bench_time.log', level=logging.INFO)


def run_bfs(datafile, expected_node_levels, alg, prt, num_dpus):

    id_str = f"{alg}_{prt}_{os.path.basename(datafile)}_{num_dpus}"
    res = f"res_{id_str}"

    make = f"NR_TASKLETS=11 BLOCK_SIZE=512 BENCHMARK_TIME=true make all"
    run = f"./bin/bfs -n {num_dpus} -a {alg} -p {prt} -o {res} {datafile}"
    subprocess.run(make, shell=True)
    process = subprocess.run(
        run, shell=True, stdout=subprocess.PIPE, encoding="utf-8")

    if process.returncode > 0:
        logging.error(f"BFS failed to complete ({id_str})")
        if os.path.exists(res):
            os.remove(res)
        return

    if not filecmp.cmp(res, expected_node_levels):
        logging.error(f"BFS output is incorrect ({id_str})")
        if os.path.exists(res):
            os.remove(res)
        return

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

    return dpu_compute_time, host_comm_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all


# dpu_compute_time, host_comm_time, pop_mram_time, fetch_res_time, total_alg, total_pop_fetch, total_all = run_bfs(
#     "data/loc-brightkite_edges.txt", "results/true/loc-brightkite_edges.txt", "src", "row", 64)

# print(f"dpu_compute_time={dpu_compute_time}, host_comm_time={host_comm_time}, pop_mram_time={pop_mram_time}, fetch_res_time={fetch_res_time}, total_alg={total_alg}, total_pop_fetch={total_pop_fetch}, total_all={total_all}")
