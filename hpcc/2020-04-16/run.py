'''
'''

import argparse, os, sys, errno, subprocess

seed_offset = 1000
num_replicates = 50

nk_config = {
    "environment_change": [
        "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 2 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 4 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 8 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 16 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 32 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 64 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 128 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 256 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 512 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 1024 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 2048 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 4096 -CHANGE_FREQUENCY 1"
        ],

    "BIT_FLIP_PROB": [
        "-BIT_FLIP_PROB 0.0001 -PHASE_2_BIT_FLIP_PROB 0.0001",
        "-BIT_FLIP_PROB 0.0003 -PHASE_2_BIT_FLIP_PROB 0.0003",
        "-BIT_FLIP_PROB 0.001 -PHASE_2_BIT_FLIP_PROB 0.001",
        "-BIT_FLIP_PROB 0.003 -PHASE_2_BIT_FLIP_PROB 0.003",
        "-BIT_FLIP_PROB 0.01 -PHASE_2_BIT_FLIP_PROB 0.01",
        "-BIT_FLIP_PROB 0.03 -PHASE_2_BIT_FLIP_PROB 0.03",
        "-BIT_FLIP_PROB 0.1 -PHASE_2_BIT_FLIP_PROB 0.1"]
}

gradient_config = {
    "environment_change": [
        "-CHANGE_MAGNITUDE 0 -CHANGE_FREQUENCY 0",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 256",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 128",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 64",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 32",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 16",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 8",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 4",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 2",
        "-CHANGE_MAGNITUDE 1 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 2 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 4 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 8 -CHANGE_FREQUENCY 1",
        "-CHANGE_MAGNITUDE 16 -CHANGE_FREQUENCY 1"
        ],

    "BIT_FLIP_PROB": [
        "-BIT_FLIP_PROB 0.0001 -PHASE_2_BIT_FLIP_PROB 0.0001",
        "-BIT_FLIP_PROB 0.0003 -PHASE_2_BIT_FLIP_PROB 0.0003",
        "-BIT_FLIP_PROB 0.001 -PHASE_2_BIT_FLIP_PROB 0.001",
        "-BIT_FLIP_PROB 0.003 -PHASE_2_BIT_FLIP_PROB 0.003",
        "-BIT_FLIP_PROB 0.01 -PHASE_2_BIT_FLIP_PROB 0.01",
        "-BIT_FLIP_PROB 0.03 -PHASE_2_BIT_FLIP_PROB 0.03",
        "-BIT_FLIP_PROB 0.1 -PHASE_2_BIT_FLIP_PROB 0.1"]
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

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--config_dir", type=str, help="Where is the configuration directory for experiment?")
    parser.add_argument("--array_id", type=int, help="Which array ID is associated with each ")

    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    array_id = args.array_id - 1 # -1 because job array starts at 1

    # Compute all combinations of NK fitness model settings and gradient fitness settings
    nk_combos = [f"{chg} {mut} -GRADIENT_MODEL 0" for chg in nk_config["environment_change"] for mut in nk_config["BIT_FLIP_PROB"] ]
    gradient_combos = [f"{chg} {mut} -GRADIENT_MODEL 1" for chg in gradient_config["environment_change"] for mut in gradient_config["BIT_FLIP_PROB"] ]

    # Combine
    combos = gradient_combos + nk_combos

    # Array ID must be valid index into combos
    if (array_id >= len(combos)):
        print("Invalid array_id,", array_id)
        exit(-1)
    condition_params = combos[array_id]
    cfg_fpath = os.path.join(config_dir, "Aagos.cfg")
    exec_fpath = os.path.join(config_dir, "Aagos")
    # Run N replicates of this condition.
    for i in range(1, num_replicates+1):
        # Compute seed for this replicate.
        seed = seed_offset + (array_id * num_replicates) + i
        # Generate run parameters, use to name run.
        run_params = condition_params + f" -SEED {seed}"
        run_name = ("RUN " + run_params).replace("-", "_").replace(" ", "_")
        print(f"Running {run_name}")
        # Make run directory
        run_dir = os.path.join(data_dir, run_name)
        mkdir_p(run_dir)
        # subprocess.run(f"cd {run_dir}", shell=True)
        # Copy configuration files into the run directory
        subprocess.run(f"cp {cfg_fpath} {run_dir}" , shell=True)
        subprocess.run(f"cp {exec_fpath} {run_dir}" , shell=True)
        # Run experiment
        run_exec = os.path.join(run_dir, "Aagos")
        run_config = os.path.join(run_dir, "Aagos.cfg")
        log_path = os.path.join(run_dir, "run.log")
        output_path = os.path.join(run_dir, "output")
        subprocess.run(f"{run_exec} {run_params} -DATA_FILEPATH {output_path} > {log_path}", shell=True)
        # Cleanup
        subprocess.run(f"rm {run_exec}" , shell=True)
        subprocess.run(f"rm {run_config}" , shell=True)
    print("Done")









if __name__ == "__main__":
    main()