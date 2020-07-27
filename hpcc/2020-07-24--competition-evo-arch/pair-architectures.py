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

run_dir_identifier = "SEED_" # All legit run directories will have this substring in their name.

update = 50000
seed_offset = 980000

run_params = "-MAX_GENS 10000 -POP_SIZE 2000 -LOAD_ANCESTOR 1"
high_mut_run_params = "-BIT_FLIP_PROB 0.1"
low_mut_run_params = "-BIT_FLIP_PROB 0.003"

job_time_request = "00:10:00"
job_memory_request = "2G"
job_name = "comp"

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
    # parser.add_argument("--src_dir", type=str, help="Where should we pull data (one or more locations)?")
    parser.add_argument("--run_dir", type=str, help="Where should we run the experiment?")
    parser.add_argument("--config_dir", type=str, help="Where should we pull data (one or more locations)?")

    # Parse command line arguments
    args = parser.parse_args()
    # src_dir = args.src_dir
    config_dir = args.config_dir
    hpcc_run_dir = args.run_dir

    dump_dir = "ancestor_pairings"

    # Make place to dump aggregated data.
    mkdir_p(dump_dir)

    # Load runs from run pairings
    run_pairings_fpath = os.path.join(config_dir, "run_pairings.csv")
    content = None
    with open(run_pairings_fpath, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    content = content[1:]
    run_pairings = [ {key:l[header_lu[key]] for key in header} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]

    # Reset these utility variables
    content = []
    header = []
    header_lu = {}

    # Extras to define:
    # - LOAD_ANCESTOR_FILE
    # - LOAD_ENV_FILE
    # (for each seed, save out ancestor ids)
    ancestor_id_lu_keys = ["seeds", "ancestor_condition", "ancestral_genome_id"]
    ancestor_id_lu = [",".join(ancestor_id_lu_keys)]
    seed_cnt = 0
    run_subs = []
    for pairing_i in range(len(run_pairings)):
        print(f"Processing pairing {pairing_i}")
        # Get 0's population (LOW MUTATION RATE)
        low_mut_src_dir = run_pairings[pairing_i]["run_dir_0"]
        low_mut_pop_path = os.path.join(low_mut_src_dir, "output", f"pop_{update}.csv")
        with open(low_mut_pop_path, "r") as fp:
            content = fp.read().strip().split("\n")
        header = content[0].split(",")
        header_lu = {header[i].strip():i for i in range(0, len(header))}
        content = content[1:]
        low_mut_pop = [ {key:l[header_lu[key]] for key in header} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
        low_mut_ancestors = [{"gene_starts":org["gene_starts"].strip("[]"), "bits":org["genome_bitstring"]} for org in low_mut_pop]

        # Get 1's population (HIGH MUTATION RATE)
        high_mut_src_dir = run_pairings[pairing_i]["run_dir_1"]
        high_mut_pop_path = os.path.join(high_mut_src_dir, "output", f"pop_{update}.csv")
        with open(high_mut_pop_path, "r") as fp:
            content = fp.read().strip().split("\n")
        header = content[0].split(",")
        header_lu = {header[i].strip():i for i in range(0, len(header))}
        content = content[1:]
        high_mut_pop = [ {key:l[header_lu[key]] for key in header} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
        high_mut_ancestors = [{"gene_starts":org["gene_starts"].strip("[]"), "bits":org["genome_bitstring"]} for org in high_mut_pop]

        # Combine ancestors
        if len(low_mut_ancestors) != len(high_mut_ancestors):
            print("# low mut ancestors =/= # high mut ancestors")
            exit(-1)
        print(f"  Number low mut ancestors: {len(low_mut_ancestors)}; Number high mut ancestors: {len(high_mut_ancestors)}")
        ancestors = []
        ancestral_genome_id = 0
        low_mut_transfer_seed = seed_offset + seed_cnt
        high_mut_transfer_seed = seed_offset + seed_cnt + 1
        seed_cnt += 2
        for i in range(len(low_mut_ancestors)):
            # Track ancestral genome id lookup info
            ancestral_id_lu_info = {
                "seeds": f"{low_mut_transfer_seed} {high_mut_transfer_seed}",
                "ancestor_condition": "LOW_MUT",
                "ancestral_genome_id": str(ancestral_genome_id)
            }
            ancestor_id_lu.append(",".join([ancestral_id_lu_info[key] for key in ancestor_id_lu_keys]))
            ancestors.append(f"{low_mut_ancestors[i]['gene_starts']},{low_mut_ancestors[i]['bits']}" )
            ancestral_genome_id += 1

        for i in range(len(high_mut_ancestors)):
            # Track ancestral genome id lookup info
            ancestral_id_lu_info = {
                "seeds": f"{low_mut_transfer_seed} {high_mut_transfer_seed}",
                "ancestor_condition": "HIGH_MUT",
                "ancestral_genome_id": str(ancestral_genome_id)
            }
            ancestor_id_lu.append(",".join([ancestral_id_lu_info[key] for key in ancestor_id_lu_keys]))
            ancestors.append(f"{high_mut_ancestors[i]['gene_starts']},{high_mut_ancestors[i]['bits']}" )
            ancestral_genome_id += 1

        # Write ancestor files out
        # - Two identical files - one for high mut and one for low mut (just to make things easier)
        low_mut_ancestor_fpath = os.path.join(config_dir, dump_dir, f"SEED_{low_mut_transfer_seed}.csv")
        high_mut_ancestor_fpath = os.path.join(config_dir, dump_dir, f"SEED_{high_mut_transfer_seed}.csv")
        with open(low_mut_ancestor_fpath, "w") as fp:
            fp.write("\n".join(ancestors))
        with open(high_mut_ancestor_fpath, "w") as fp:
            fp.write("\n".join(ancestors))

        # Compose run information for low and high mutation rate treatments
        load_env_file = run_pairings[pairing_i]["env"]
        gradient_model = run_pairings[pairing_i]["gradient_model"]

        # -- low mut --
        low_mut_params = f"{run_params} {low_mut_run_params} -SEED {low_mut_transfer_seed} -LOAD_ENV_FILE {load_env_file} -LOAD_ANCESTOR_FILE {low_mut_ancestor_fpath} -GRADIENT_MODEL {gradient_model}"
        low_mut_run_dir = os.path.join(hpcc_run_dir, f"SEED_{low_mut_transfer_seed}")
        run_subs.append({
            "run_params": low_mut_params,
            "run_dir": low_mut_run_dir
        })

        # -- high mut --
        high_mut_params = f"{run_params} {high_mut_run_params} -SEED {high_mut_transfer_seed} -LOAD_ENV_FILE {load_env_file} -LOAD_ANCESTOR_FILE {high_mut_ancestor_fpath} -GRADIENT_MODEL {gradient_model}"
        high_mut_run_dir = os.path.join(hpcc_run_dir, f"SEED_{high_mut_transfer_seed}")
        run_subs.append({
            "run_params": high_mut_params,
            "run_dir": high_mut_run_dir
        })


    # Write out lookup table
    with open("ancestor_id_lu.csv", "w") as fp:
        fp.write("\n".join(ancestor_id_lu))

    print("Runs to submit: " + str(len(run_subs)))
    resub_logic = ""
    array_id = 1
    for resub in run_subs:
        sub_params = resub["run_params"]
        run_logic = base_run_logic
        run_logic = run_logic.replace("<<RESUB_ID>>", str(array_id))
        run_logic = run_logic.replace("<<RUN_DIR>>", resub["run_dir"])
        run_logic = run_logic.replace("<<RUN_PARAMS>>", f"'{sub_params}'")

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