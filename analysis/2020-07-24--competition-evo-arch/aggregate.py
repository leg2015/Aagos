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
import argparse, os, sys, errno, csv

run_dir_identifier = "SEED_" # All legit run directories will have this substring in their name.

config_exclude = {"LOAD_ANCESTOR_FILE", "PHASE_2_ENV_FILE", "DATA_FILEPATH", "SNAPSHOT_INTERVAL", "PRINT_INTERVAL", "SUMMARY_INTERVAL"}
gene_stats_exclude = {}
rep_org_exclude = {f"site_cnt_{i}_gene_occupancy" for i in range(0, 128)}
rep_org_exclude.add("gene_neighbors")

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
    parser.add_argument("--data", type=str, nargs="+", help="Where should we pull data (one or more locations)?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--update", type=int, help="Which update should we pull?")
    parser.add_argument("--ancestor_map", type=str,help="What ancestor ids belong with each original environment?")

    # Parse command line arguments
    args = parser.parse_args()
    data_dirs = args.data
    dump_dir = args.dump
    update = args.update
    ancestor_map_fpath = args.ancestor_map


    # Load ancestor map
    # - For each run seed, map each genome id to high mut or low mut
    #   {run-seed: {genome_id: 'high', ...}}
    content = None
    with open(ancestor_map_fpath, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    header_lu = {header[i].strip():i for i in range(0, len(header))}
    content = content[1:]
    ancestor_map_content = [ {key:l[header_lu[key]] for key in header} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]

    ancestor_genome_lu = {}
    for line in ancestor_map_content:
        seeds = line["seeds"].split(" ")
        for seed in seeds:
            if not seed in ancestor_genome_lu: ancestor_genome_lu[seed] = {}
            ancestor_genome_lu[seed][line["ancestral_genome_id"]] = line["ancestor_condition"]

    # Load data...

    # Are all data directories for real?
    if any([not os.path.exists(loc) for loc in data_dirs]):
        print("Unable to locate all data directories. Locatability:", {loc: os.path.exists(loc) for loc in data_dirs})
        exit(-1)
    # Make place to dump aggregated data.
    mkdir_p(dump_dir)
    # Aggregate a list of all runs.
    run_dirs = [os.path.join(data_dir, run_dir) for data_dir in data_dirs for run_dir in os.listdir(data_dir) if run_dir_identifier in run_dir]
    # Sort run directories by seed
    run_dirs.sort(key=lambda x : int(x.split("_")[-1]))
    print(f"Found {len(run_dirs)} run directories.")

    combined_header_set = set()
    aggregate_info = []

    # aggregate info:
    # - # high mut lineages
    # - # low mut lineages
    # - dominant_lineage_type
    # - run settings

    for run in run_dirs:
        print(f"Extracting information from {run}")
        run_config_path = os.path.join(run, "output", "run_config.csv")
        pop_snapshot_path = os.path.join(run, "output", f"pop_{update}.csv")

        # Extract the run settings
        if not os.path.exists(run_config_path):
            print(f"Failed to open run configuration file: {run_config_path}")
            exit(-1)
        run_settings = extract_settings(run_config_path)

        run_seed = run_settings["SEED"]

        # Load pop snapshot
        snapshot_content = None
        with open(pop_snapshot_path) as fp:
            snapshot_content = fp.read().strip().split("\n")
        snapshot_header = snapshot_content[0].split(",")
        snapshot_content = snapshot_content[1:]
        pop = [{snapshot_header[i]: l[i] for i in range(0, len(l))} for l in csv.reader(snapshot_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]

        low_mut_arch_cnt = 0
        high_mut_arch_cnt = 0

        for org in pop:
            org_arch = ancestor_genome_lu[run_seed][org["ancestral_id"]]
            low_mut_arch_cnt += int(org_arch == "LOW_MUT")
            high_mut_arch_cnt += int(org_arch == "HIGH_MUT")

        print(f"High mut arch cnt: {high_mut_arch_cnt}")
        print(f"Low mut arch cnt: {low_mut_arch_cnt}")

        dom_lineage_type = "NONE"
        if low_mut_arch_cnt > high_mut_arch_cnt: dom_lineage_type = "LOW_MUT"
        if high_mut_arch_cnt > low_mut_arch_cnt: dom_lineage_type = "HIGH_MUT"

        extra_field_keys = ["high_mut_lineage_cnt", "low_mut_lineage_cnt", "dominant_lineage_type"]
        extra_fields = {
            "high_mut_lineage_cnt": high_mut_arch_cnt,
            "low_mut_lineage_cnt": low_mut_arch_cnt,
            "dominant_lineage_type": dom_lineage_type
        }

        config_fields = [field for field in run_settings if (field not in config_exclude)]
        fields = [f for f in extra_field_keys] + config_fields
        combined_header_set.add(",".join(fields))
        if len(combined_header_set) > 1:
            print("Header mismatch!")
            exit(-1)
        line = [extra_fields[f] for f in extra_field_keys] + [run_settings[s] for s in config_fields]
        aggregate_info.append(line)

    # Output aggregate information.
    out_content = list(combined_header_set)[0] + "\n"
    out_content += "\n".join([",".join(map(str, line)) for line in aggregate_info])
    out_path = os.path.join(dump_dir, "agg_data.csv")
    with open(out_path, "w") as fp:
        fp.write(out_content)
    print(f"Done! Output written to {out_path}")








if __name__ == "__main__":
    main()