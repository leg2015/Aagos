'''
Search data directory for incomplete runs, generate submission script for all
incomplete (unfinished or not started) runs.
'''

import argparse, os, sys, errno, subprocess, csv

seed_offset = 60000
default_num_replicates = 100

shared_config = {

    "environment_change": [
        "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0",
    ],

    "BIT_FLIP_PROB": [
        "-BIT_FLIP_PROB 0.0001 -PHASE_2_BIT_FLIP_PROB 0.0001",
        "-BIT_FLIP_PROB 0.0003 -PHASE_2_BIT_FLIP_PROB 0.0003",
        "-BIT_FLIP_PROB 0.001 -PHASE_2_BIT_FLIP_PROB 0.001",
        "-BIT_FLIP_PROB 0.003 -PHASE_2_BIT_FLIP_PROB 0.003",
        "-BIT_FLIP_PROB 0.01 -PHASE_2_BIT_FLIP_PROB 0.01",
        "-BIT_FLIP_PROB 0.03 -PHASE_2_BIT_FLIP_PROB 0.03",
        "-BIT_FLIP_PROB 0.1 -PHASE_2_BIT_FLIP_PROB 0.1"
    ],

    "GENOME_SIZE": [
        "-NUM_BITS 128 -NUM_GENES 16 -MAX_SIZE 1024",
        "-NUM_BITS 64 -NUM_GENES 8 -MAX_SIZE 512",
        "-NUM_BITS 256 -NUM_GENES 32 -MAX_SIZE 2048"
    ],

    "SELECTION": [
        "-TOURNAMENT_SIZE 1 -PHASE_2_TOURNAMENT_SIZE 8"
    ]

}



'''
This is functionally equivalent to the mkdir -p [fname] bash command
'''
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def extract_settings(run_config_path):
    content = None
    with open(run_config_path, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    content = content[1:]
    configs = [l for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
    return {param[header_lu["parameter"]]:param[header_lu["value"]] for param in configs}


def is_run_complete(path):
    # (1) Does the run directory exist?
    if not os.path.exists(path): return False
    # (2) If the run directory exists, did the run complete?
    #     Is there a run config file?
    run_config_path = os.path.join(path, "output", "run_config.csv")
    if not os.path.exists(run_config_path): return False
    #    The run config file exists, extract parameters.
    run_params = extract_settings(run_config_path)
    final_gen = run_params["TOTAL_GENS"] # We'll look for this generation in the fitness.csv file
    fitness_file_path = os.path.join(path, "output", "fitness.csv")
    if not os.path.exists(fitness_file_path): return False
    fitness_contents = None
    with open(fitness_file_path, "r") as fp:
        fitness_contents = fp.read().strip().split("\n")
    if len(fitness_contents) == 0: return False
    header = fitness_contents[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    if len(header_lu) != fitness_contents[-1].split(","): return False
    if fitness_contents[-1].split(",")[header_lu["update"]] != final_gen: return False
    return True

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--config_dir", type=str, help="Where is the configuration directory for experiment?")
    # parser.add_argument("--array_id", type=int, help="Which array ID is associated with each ")
    parser.add_argument("--replicates", type=int, default=default_num_replicates, help="How many replicates should we run of each condition?")
    parser.add_argument("--query_condition_cnt", action="store_true", help="How many conditions?")

    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    # array_id = args.array_id - 1        # -1 because job array starts at 1
    num_replicates = args.replicates

    # Compute all combinations of NK fitness model settings and gradient fitness settings
    combos = [f"{chg} {mut} {genome} {fitness}" for fitness in ["-GRADIENT_MODEL 0", "-GRADIENT_MODEL 1"] for chg in shared_config["environment_change"] for mut in shared_config["BIT_FLIP_PROB"] for genome in shared_config["GENOME_SIZE"] ]
    if (args.query_condition_cnt):
        print("Conditions", combos)
        print(f"Number of conditions: {len(combos)}")
        exit(0)

    cfg_fpath = os.path.join(config_dir, "Aagos.cfg")
    exec_fpath = os.path.join(config_dir, "Aagos")

    # Find complete/incomplete runs.
    num_finished = 0
    resubmissions = []
    for condition_id in range(0, len(combos)):
        condition_params = combos[condition_id]
        print(f"Processing condition: {condition_params}")
        # Run N replicates of this condition.
        for i in range(1, num_replicates+1):
            # Compute seed for this replicate.
            seed = seed_offset + (condition_id * num_replicates) + i
            # Generate run parameters, use to name run.
            run_params = condition_params + f" -SEED {seed}"
            run_name = ("RUN " + run_params).replace("-", "_").replace(" ", "_")
            run_dir = os.path.join(data_dir, run_name)
            # (1) Does the run directory exist?
            run_complete = is_run_complete(run_dir)
            num_finished += int(run_complete)
            # if not run_complete: resubmissions.append({"run_dir": run_dir, "run_params": run_params})
            if not run_complete: resubmissions.append(run_params)

            # mkdir_p(run_dir)
            # subprocess.run(f"cd {run_dir}", shell=True)
            # Copy configuration files into the run directory
            # subprocess.run(f"cp {cfg_fpath} {run_dir}" , shell=True)
            # subprocess.run(f"cp {exec_fpath} {run_dir}" , shell=True)
            # Run experiment
            # run_exec = os.path.join(run_dir, "Aagos")
            # run_config = os.path.join(run_dir, "Aagos.cfg")
            # log_path = os.path.join(run_dir, "run.log")
            # output_path = os.path.join(run_dir, "output")
            # subprocess.run(f"{run_exec} {run_params} -DATA_FILEPATH {output_path} > {log_path}", shell=True)
            # Cleanup
            # subprocess.run(f"rm {run_exec}" , shell=True)
            # subprocess.run(f"rm {run_config}" , shell=True)
    print(f"Runs finished: {num_finished}")
    print(f"Resubmissions: {len(resubmissions)}")
    for resub in resubmissions:
        print(f"  - {resub}")

if __name__ == "__main__":
    main()