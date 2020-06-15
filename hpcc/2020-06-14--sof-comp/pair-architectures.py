'''
Aggregate data from each Aagos run.
For each run, collect information from both phase 1 and 2 of evolution (last update for each)

Gene stats info (prepend 'pop_'):
- update,evo_phase,
- mean_neutral_sites,min_neutral_sites,max_neutral_sites,variance_neutral_sites,mean_single_gene_sites,min_single_gene_sites,max_single_gene_sites,variance_single_gene_sites,mean_multi_gene_sites,min_multi_gene_sites,max_multi_gene_sites,variance_multi_gene_sites,mean_site_occupancy,min_site_occupancy,max_site_occupancy,variance_site_occupancy,mean_neighbor_genes,min_neighbor_genes,max_neighbor_genes,variance_neighbor_genes,mean_coding_sites,min_coding_sites,max_coding_sites,variance_coding_sites,mean_genome_length,min_genome_length,max_genome_length,variance_genome_length

Representative info (prepend 'rep_'):
- update,evo_phase,
- fitness,genome_length,coding_sites,neutral_sites,gene_starts,genome_bitstring,gene_size,gene_neighbors,avg_gene_neighbors,site_cnt_0_gene_occupancy,site_cnt_1_gene_occupancy,site_cnt_2_gene_occupancy,site_cnt_3_gene_occupancy,site_cnt_4_gene_occupancy,site_cnt_5_gene_occupancy,site_cnt_6_gene_occupancy,site_cnt_7_gene_occupancy,site_cnt_8_gene_occupancy,site_cnt_9_gene_occupancy,site_cnt_10_gene_occupancy,site_cnt_11_gene_occupancy,site_cnt_12_gene_occupancy,site_cnt_13_gene_occupancy,site_cnt_14_gene_occupancy,site_cnt_15_gene_occupancy,site_cnt_16_gene_occupancy

Config info: [collect all]
'''
import argparse, os, sys, errno, csv, copy

run_dir_identifier = "__SEED_" # All legit run directories will have this substring in their name.

runs_per_pair = 100
seed_offset = 900000
update = 50000
# update = 200
job_time_request = "00:10:00"
job_memory_request = "2G"
job_name = "cmp"

categories = {
  "LMST-GRD":{"GRADIENT_MODEL":"1", "BIT_FLIP_PROB":"0.003", "CHANGE_MAGNITUDE":"0", "CHANGE_FREQUENCY":"0"},
  "HMST-GRD":{"GRADIENT_MODEL":"1", "BIT_FLIP_PROB":"0.1", "CHANGE_MAGNITUDE":"0", "CHANGE_FREQUENCY":"0"},
  "LMCH-GRD":{"GRADIENT_MODEL":"1", "BIT_FLIP_PROB":"0.003", "CHANGE_MAGNITUDE":"1", "CHANGE_FREQUENCY":"4"},
  "HMCH-GRD":{"GRADIENT_MODEL":"1", "BIT_FLIP_PROB":"0.1", "CHANGE_MAGNITUDE":"1", "CHANGE_FREQUENCY":"4"},
  "LMST-NK":{"GRADIENT_MODEL":"0", "BIT_FLIP_PROB":"0.003", "CHANGE_MAGNITUDE":"0", "CHANGE_FREQUENCY":"0"},
  "HMST-NK":{"GRADIENT_MODEL":"0", "BIT_FLIP_PROB":"0.1", "CHANGE_MAGNITUDE":"0", "CHANGE_FREQUENCY":"0"},
  "LMCH-NK":{"GRADIENT_MODEL":"0", "BIT_FLIP_PROB":"0.003", "CHANGE_MAGNITUDE":"64", "CHANGE_FREQUENCY":"1"},
  "HMCH-NK":{"GRADIENT_MODEL":"0", "BIT_FLIP_PROB":"0.1", "CHANGE_MAGNITUDE":"64", "CHANGE_FREQUENCY":"1"}
}

# NOTE: First architecture in file always from low mutation rate environment
pairings = [
    ["LMST-GRD", "HMST-GRD"],
    ["LMCH-GRD", "HMCH-GRD"],
    ["LMCH-GRD", "HMST-GRD"],

    ["LMST-NK", "HMST-NK"],
    ["LMCH-NK", "HMCH-NK"],
    ["LMCH-NK", "HMST-NK"]
]

high_mut = "-BIT_FLIP_PROB 0.1"
low_mut = "-BIT_FLIP_PROB 0.003"

base_resub_script = \
"""#!/bin/bash
########## Define Resources Needed with SBATCH Lines ##########

#SBATCH --time=<<TIME_REQUEST>>          # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --array=<<ARRAY_ID_RANGE>>
#SBATCH --mem=<<MEMORY_REQUEST>>        # memory required per node - amount of memory (in bytes)
#SBATCH --job-name <<JOB_NAME>>         # you can give your job a name for easier identification (same as -J)
#SBATCH --account=devolab

########## Command Lines to Run ##########

EXEC=Aagos
CONFIG_DIR=<<CONFIG_DIR>>

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.7.0

<<RESUBMISSION_LOGIC>>

mkdir -p ${RUN_DIR}
cd ${RUN_DIR}
cp ${CONFIG_DIR}/Aagos.cfg .
cp ${CONFIG_DIR}/${EXEC} .

./${EXEC} ${RUN_PARAMS} > run.log

rm Aagos.cfg
rm ${EXEC}

"""

base_run_logic = \
"""
if [[ ${SLURM_ARRAY_TASK_ID} -eq <<RESUB_ID>> ]] ; then
    RUN_DIR=<<RUN_DIR>>
    RUN_PARAMS=<<RUN_PARAMS>>
fi
"""


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

'''
Given the path to a run's config file, extract the run's settings.
'''
def extract_settings(run_config_path):
    content = None
    with open(run_config_path, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    content = content[1:]
    configs = [l for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
    return {param[header_lu["parameter"]]:param[header_lu["value"]] for param in configs}

def main():
    # Setup the command line argument parser
    parser = argparse.ArgumentParser(description="Data aggregation script")
    parser.add_argument("--data_dir", type=str, help="Where should we pull data (one or more locations)?")
    parser.add_argument("--run_dir", type=str, help="Where should we run the experiment?")
    parser.add_argument("--config_dir", type=str, help="Where should we pull data (one or more locations)?")

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    hpcc_run_dir = args.run_dir
    dump_dir = "ancestor_pairings"

    # Make place to dump aggregated data.
    mkdir_p(data_dir)
    mkdir_p(dump_dir)

    # Aggregate a list of all runs.
    run_dirs = [os.path.join(data_dir, run_dir) for run_dir in os.listdir(data_dir) if run_dir_identifier in run_dir]

    # Sort run directories by seed
    run_dirs.sort(key=lambda x : int(x.split("_")[-1]))
    print(f"Found {len(run_dirs)} run directories.")

    run_orgs = [None for _ in run_dirs]
    run_configurations = [None for _ in run_dirs]
    run_ids_by_category = {cat:[] for cat in categories}

    # (1) Categorize each run
    for run_index in range(0, len(run_dirs)):
        run_dir = run_dirs[run_index]
        print("===")
        print(f"Categorizing {run_dir}")
        run_config_path = os.path.join(run_dir, "output", "run_config.csv")
        rep_org_path = os.path.join(run_dir, "output", "representative_org.csv")

        # Extract the run settings
        if not os.path.exists(run_config_path):
            print(f"Failed to open run configuration file: {run_config_path}")
            exit(-1)
        run_settings = extract_settings(run_config_path)

        # Extract the representative organism content
        content = None
        with open(rep_org_path) as fp:
            content = fp.read().strip().split("\n")
        rep_org_header = content[0].split(",")
        rep_org_header_lu = {rep_org_header[i].strip():i for i in range(0, len(rep_org_header))}
        content = content[1:]
        rep_org = {int(l[rep_org_header_lu["update"]]):{rep_org_header[i]: l[i] for i in range(0, len(l))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) if int(l[rep_org_header_lu["update"]]) == update}
        print("There should only be 1 rep org..", len(rep_org))
        rep_org = rep_org[update]
        # if "gene_starts" in rep_org:
        #     rep_org["gene_starts"] = "\"" + rep_org["gene_starts"] + "\""

        # Store organism info for this run
        run_orgs[run_index] = copy.deepcopy(rep_org)
        run_configurations[run_index] = copy.deepcopy(run_settings)

        # Figure out which category this run falls into
        run_category = None
        for category_name in categories:
            category_params = categories[category_name]
            if all([category_params[param]==run_settings[param] for param in category_params]):
                run_category = category_name
                break
        run_ids_by_category[run_category].append(run_index)
    print("Run ids by category:" + str(run_ids_by_category))

    # (2) Form ancestor pairings
    #     - For each pairing type, generate N pairs
    #     - For each of N pairs, generate a
    run_subs = []
    cnt = 0
    for pairing_type in pairings:
        cat_a = pairing_type[0]
        cat_b = pairing_type[1]
        pairing_dump = os.path.join(dump_dir, "_".join(pairing_type))
        mkdir_p(pairing_dump)
        # Make a directory for this pairing
        if len(run_ids_by_category[cat_a]) < runs_per_pair:
            print(f"Insufficient runs for category {cat_a}")
            exit(-1)
        if len(run_ids_by_category[cat_b]) < runs_per_pair:
            print(f"Insufficient runs for category {cat_b}")
            exit(-1)
        cat_a_run_ids = []
        cat_b_run_ids = []
        while len(cat_a_run_ids) < runs_per_pair:
            cat_a_run_ids.append(run_ids_by_category[cat_a].pop())
        while len(cat_b_run_ids) < runs_per_pair:
            cat_b_run_ids.append(run_ids_by_category[cat_b].pop())
        for cat_a_run_id, cat_b_run_id in zip(cat_a_run_ids, cat_b_run_ids):
            # Extract architectures
            cat_a_org = run_orgs[cat_a_run_id]
            cat_b_org = run_orgs[cat_b_run_id]
            cat_a_settings = run_configurations[cat_a_run_id]
            cat_b_settings = run_configurations[cat_b_run_id]
            if (cat_a_settings["GRADIENT_MODEL"] != cat_b_settings["GRADIENT_MODEL"]):
                print("Fitness model mismatch!")
                exit(-1)
            cat_a_bits = cat_a_org["genome_bitstring"].replace("1", "0")
            cat_b_bits = cat_b_org["genome_bitstring"].replace("1", "0")
            arch_content = cat_a_org["gene_starts"].strip("[]") + "," + cat_a_bits + "\n"
            arch_content += cat_b_org["gene_starts"].strip("[]") + "," + cat_b_bits
            # Write architectures to a file
            arch_fname = "_".join(pairing_type) + "_" + str(cat_a_run_id) + "-" + str(cat_b_run_id) + ".csv"
            arch_fpath = os.path.join(pairing_dump, arch_fname)
            with open(arch_fpath, "w") as fp:
                fp.write(arch_content)
            # Generate a stubs for this run (append to run_subs)
            #  - run stub should have command line configuration info
            #    - low mutation rate
            #    - high mutation
            params = {}
            params["GRADIENT_MODEL"] = cat_a_settings["GRADIENT_MODEL"]
            params["CHANGE_MAGNITUDE"] = "0"
            params["CHANGE_FREQUENCY"] = "0"
            params["MAX_GENS"] = "1000"
            params["LOAD_ANCESTOR"] = "1"
            params["LOAD_ANCESTOR_FILE"] = os.path.join(config_dir, arch_fpath)
            params["GENE_MOVE_PROB"] = "0.0"
            params["BIT_INS_PROB"] = "0.0"
            params["BIT_DEL_PROB"] = "0.0"
            params["SUMMARY_INTERVAL"] = "50"
            params["SNAPSHOT_INTERVAL"] = "1000"
            params["SEED"] = seed_offset + int(cnt)
            cnt += 2
            # Low mutation rate competition
            low_mut_stub = " ".join([f"-{cfg} {params[cfg]}" for cfg in params] + [low_mut])
            run_name = "_".join(pairing_type) + "__LOW_MUT__" + "SEED_" + str(params["SEED"])
            run_subs.append({
                "run_params": low_mut_stub,
                "run_dir": os.path.join(hpcc_run_dir, run_name)
            })
            # High mutation rate competition
            params["SEED"] += 1
            high_mut_stub = " ".join([f"-{cfg} {params[cfg]}" for cfg in params] + [high_mut])
            run_name = "_".join(pairing_type) + "__HIGH_MUT__" + "SEED_" + str(params["SEED"])
            run_subs.append({
                "run_params": high_mut_stub,
                "run_dir": os.path.join(hpcc_run_dir, run_name)
            })

    # (3) Generate submission script from run subs
    print("Runs to submit: " + str(len(run_subs)))
    resub_logic = ""
    array_id = 1
    for resub in run_subs:
        run_params = resub["run_params"]
        run_logic = base_run_logic
        run_logic = run_logic.replace("<<RESUB_ID>>", str(array_id))
        run_logic = run_logic.replace("<<RUN_DIR>>", resub["run_dir"])
        run_logic = run_logic.replace("<<RUN_PARAMS>>", f"'{run_params}'")

        resub_logic += run_logic
        array_id += 1

    script = base_resub_script
    script = script.replace("<<TIME_REQUEST>>", job_time_request)
    script = script.replace("<<ARRAY_ID_RANGE>>", f"1-{len(run_subs)}")
    script = script.replace("<<MEMORY_REQUEST>>", job_memory_request)
    script = script.replace("<<JOB_NAME>>", job_name)
    script = script.replace("<<CONFIG_DIR>>", config_dir)
    script = script.replace("<<RESUBMISSION_LOGIC>>", resub_logic)

    with open("competition_sub.sb", "w") as fp:
        fp.write(script)


if __name__ == "__main__":
    main()