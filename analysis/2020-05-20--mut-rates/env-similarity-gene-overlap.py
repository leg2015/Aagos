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
import argparse, os, sys, errno, csv, itertools

run_dir_identifier = "__SEED_" # All legit run directories will have this substring in their name.

config_exclude = {"PHASE_2_ENV_FILE", "DATA_FILEPATH", "SNAPSHOT_INTERVAL", "PRINT_INTERVAL", "SUMMARY_INTERVAL"}
gene_stats_exclude = {}
rep_org_exclude = {"gene_starts", "genome_bitstring", "gene_neighbors"}

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

'''
Find the max (hamming) similarity between a and b for all possible shifted alignments

E.g.,

a = 0123
b = 4567

   0123
4567
 4567
  4567
   4567
    4567
     4567
      4567
'''
def max_aligned_similarity(a, b):
    if len(a) != len(b):
        print("Must be equal length!")
        exit(-1)
    # perfectly aligned
    num_bits = len(a)
    max_similarity = -1

    for shift in range(0, num_bits):
        right_alignment_similarity = 0
        left_alignment_similarity = 0

        for i in range(0, len(a)-shift):
            right_alignment_similarity += int(a[i+shift] == b[i])

        if shift == 0:
            max_similarity = right_alignment_similarity
            continue

        for i in range(0, len(b)-shift):
            left_alignment_similarity += int(a[i] == b[i+shift])

        max_similarity = max(max_similarity, right_alignment_similarity, left_alignment_similarity)

    return max_similarity

def main():
    # Setup the command line argument parser
    parser = argparse.ArgumentParser(description="Data aggregation script")
    parser.add_argument("--data", type=str, nargs="+", help="Where should we pull data (one or more locations)?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--updates", type=int, nargs="+", help="Which updates should we pull?")
    # Parse command line arguments
    args = parser.parse_args()
    data_dirs = args.data
    dump_dir = args.dump
    updates = args.updates

    # Are all data directories for real?
    if any([not os.path.exists(loc) for loc in data_dirs]):
        print("Unable to locate all data directories. Locatability:", {loc: os.path.exists(loc) for loc in data_dirs})
        exit(-1)
    # Is there at least one update in updates?
    # updates = set(updates)
    if len(updates) == 0:
        print("No target updates provided.")
        exit(-1)

    # Make place to dump aggregated data.
    mkdir_p(dump_dir)
    # Aggregate a list of all runs.
    run_dirs = [os.path.join(data_dir, run_dir) for data_dir in data_dirs for run_dir in os.listdir(data_dir) if run_dir_identifier in run_dir]
    # Sort run directories by seed
    run_dirs.sort(key=lambda x : int(x.split("_")[-1]))
    print(f"Found {len(run_dirs)} run directories.")

    rep_org_header_set = set() # We'll use this as a sanity check for guaranteeing that all file headers match.
    # gene_stats_header_set = set()
    combined_header_set = set()
    # aggregate_header = []
    aggregate_info = []
    for run in run_dirs:
        print(f"Extracting information from {run}")
        run_config_path = os.path.join(run, "output", "run_config.csv")
        rep_org_path = os.path.join(run, "output", "representative_org.csv")
        # gene_stats_path = os.path.join(run, "output", "gene_stats.csv")
        env_path = os.path.join(run, "output", "environment.csv")
        # Extract the run settings
        if not os.path.exists(run_config_path):
            print(f"Failed to open run configuration file: {run_config_path}")
            exit(-1)
        run_settings = extract_settings(run_config_path)

        # Only aggregate from static environment, gradient model runs
        if run_settings["GRADIENT_MODEL"] != "1": continue

        # Extract gene stats content
        # content = None
        # with open(gene_stats_path, "r") as fp:
        #     content = fp.read().strip().split("\n")
        # gene_stats_header = content[0].split(",")
        # gene_stats_header_lu = {gene_stats_header[i].strip():i for i in range(0, len(gene_stats_header))}
        # content = content[1:]

        # Collect the gene stats that matter
        # gene_stats = {int(l[gene_stats_header_lu["update"]]):{gene_stats_header[i]: l[i] for i in range(0, len(l))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) if int(l[gene_stats_header_lu["update"]]) in updates}

        # Check that this gene stats header matches previous.
        # gene_stats_header_set.add(",".join(gene_stats_header))
        # if len(gene_stats_header_set) > 1:
        #     print(f"Header mismatch! ({gene_stats_path})")
        #     exit(-1)

        # Extract representative organism content
        content = None
        with open(rep_org_path, "r") as fp:
            content = fp.read().strip().split("\n")
        rep_org_header = content[0].split(",")
        rep_org_header_lu = {rep_org_header[i].strip():i for i in range(0, len(rep_org_header))}
        content = content[1:]
        rep_org = {int(l[rep_org_header_lu["update"]]):{rep_org_header[i]: l[i] for i in range(0, len(l))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) if int(l[rep_org_header_lu["update"]]) in updates}

        # Extract environment content
        env_content = None
        with open(env_path, "r") as fp:
            env_content = fp.read().strip().split("\n")
        env_header = env_content[0].split(",")
        env_header_lu = {env_header[i].strip():i for i in range(0, len(env_header))}
        env_content = env_content[1:]
        envs = {int(l[env_header_lu["update"]]): {env_header[i] : l[i] for i in range(0, len(l))} for l in csv.reader(env_content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) if int(l[env_header_lu["update"]]) in updates}
        # Measure pairwise similarity between each gene target
        all_pairs = [tuple(p) for p in itertools.combinations(range(0, int(run_settings["NUM_GENES"])), 2)]

        if len(envs) != len(updates):
            print(f"Failed to find all environments for all updates ({str(updates)})")
            exit(-1)

        # Check that this representative org header matches previous.
        rep_org_header_set.add(",".join(rep_org_header))
        if len(rep_org_header_set) > 1:
            print(f"Header mismatch! ({rep_org_header})")
            exit(-1)


        # Build a joint header.
        # - gene stats fields
        # gene_stats_fields = [field for field in gene_stats_header if field not in gene_stats_exclude]
        # field_set = set(gene_stats_fields)
        field_set = set()
        # - rep org fields
        rep_org_fields = [field for field in rep_org_header if (field not in rep_org_exclude) and (field not in field_set)]
        field_set.union(set(rep_org_fields))
        # - run configuration fields
        config_fields = [field for field in run_settings if (field not in config_exclude) and (field not in field_set)]
        fields = rep_org_fields + config_fields + ["gene_pair", "gene_a", "gene_b", "gene_pair_overlap", "gene_pair_target_similarity", "gene_pair_target_max_alignment_similarity"]
        # Combine all fields into single header (double check that this matches all previous computed headers)
        combined_header_set.add(",".join(fields))
        if len(combined_header_set) > 1:
            print("Header mismatch!")
            exit(-1)

        run_config_line_comp = [run_settings[field] for field in config_fields]
        for update in updates:
            org = rep_org[update]
            env = envs[update]
            env_targets = env["env_state"].strip("[]").strip().split(" ")
            # Calculate pairwise environment similarity
            env_similarity = {pair:sum( [int(a == b) for a,b in zip(env_targets[pair[0]], env_targets[pair[1]]) ]) for pair in all_pairs}
            env_max_alignment_similarity = {pair: max_aligned_similarity( env_targets[pair[0]], env_targets[pair[1]] ) for pair in all_pairs }
            # env_similarity_normalized = {pair: env_similarity[pair] / float(run_settings["GENE_SIZE"]) for pair in all_pairs}

            # for pair in all_pairs:
            #     a = env_targets[pair[0]]
            #     b = env_targets[pair[1]]
            #     print(f"Pair: {str(pair)}")
            #     print(f"  a: {a}")
            #     print(f"  b: {b}")
            #     print(f"  similarity = {env_similarity[pair]}")

            # Collect pairwise gene overlap
            # - for each gene, compute sites => set, compute union
            gene_starts = list(map(int, org["gene_starts"].strip("[]").split(",")))
            genome_len = int(org["genome_length"])
            gene_size = int(run_settings["GENE_SIZE"])
            gene_sites = [set([(gene_starts[gene_id] + i) % genome_len for i in range(0, gene_size)]) for gene_id in range(len(gene_starts))]
            gene_overlap = {pair: len(gene_sites[pair[0]] & gene_sites[pair[1]]) for pair in all_pairs}
            # print(gene_starts)
            # print(gene_sites)
            # print(gene_overlap)
            for pair in all_pairs:
                # ["gene_pair", "gene_a", "gene_b", "gene_pair_overlap", "gene_pair_target_similarity"]
                line = [rep_org[update][field] for field in rep_org_fields] + run_config_line_comp + [f"\"({pair[0]} {pair[1]})\"", str(pair[0]), str(pair[1]), str(gene_overlap[pair]), str(env_similarity[pair]), str(env_max_alignment_similarity[pair])]
                aggregate_info.append(line)



            # Grab (gene_stats line content for this update) + (rep_org line content for this update) + (run configuration for this update)


    # Output aggregate information.
    print(combined_header_set)
    out_content = list(combined_header_set)[0] + "\n"
    out_content += "\n".join([",".join(map(str, line)) for line in aggregate_info])
    out_path = os.path.join(dump_dir, "env_sim_gene_overlap.csv")
    with open(out_path, "w") as fp:
        fp.write(out_content)
    print(f"Done! Output written to {out_path}")

if __name__ == "__main__":
    main()