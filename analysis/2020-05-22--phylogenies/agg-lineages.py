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

run_dir_identifier = "__SEED_" # All legit run directories will have this substring in their name.

config_exclude = {"LOAD_ANCESTOR_FILE", "PHASE_2_ENV_FILE", "DATA_FILEPATH", "SNAPSHOT_INTERVAL", "PRINT_INTERVAL", "SUMMARY_INTERVAL"}

rep_org_exclude = {"gene_neighbors"}



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
    # Parse command line arguments
    args = parser.parse_args()
    data_dirs = args.data
    dump_dir = args.dump
    update = args.update

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

    lineage_summary_head_set = set()
    lineage_summary_info = []
    full_lineage_info = [] # WARNING - this will get HUGE
    for run in run_dirs:
        print(f"Extracting information from {run}")
        run_config_path = os.path.join(run, "output", "run_config.csv")

        rep_org_path = os.path.join(run, "output", "representative_org.csv")
        phylogeny_path = os.path.join(run, "output", f"phylo_{update}.csv")

        # Extract the run settings
        if not os.path.exists(run_config_path):
            print(f"Failed to open run configuration file: {run_config_path}")
            exit(-1)
        run_settings = extract_settings(run_config_path)

        # Extract representative organism content
        content = None
        with open(rep_org_path, "r") as fp:
            content = fp.read().strip().split("\n")

        rep_org_header = content[0].split(",")
        rep_org_header_lu = {rep_org_header[i].strip():i for i in range(0, len(rep_org_header))}
        content = content[1:]
        rep_org = {int(l[rep_org_header_lu["update"]]):{rep_org_header[i]: l[i] for i in range(0, len(l))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) if int(l[rep_org_header_lu["update"]]) == update}
        rep_org = rep_org[update]


        # Pull out the lineage
        content = None
        with open(phylogeny_path, "r") as fp:
            content = fp.read().strip().split("\n")
        phylo_header = content[0].split(",")
        phylo_header_lu = {phylo_header[i].strip():i for i in range(len(phylo_header))}
        content = content[1:]
        phylogeny = { l[phylo_header_lu["id"]]: {phylo_header[i]: l[i] for i in range(len(l))}  for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) }

        # Find representative organism in phylogeny (to start lineage from)
        tail = None
        cur_tail_origin = 0
        for ancestor_id in phylogeny:
            fossil = phylogeny[ancestor_id]
            is_candidate = (fossil["genome_bitstring"] == rep_org["genome_bitstring"]) and (fossil["gene_starts"] == rep_org["gene_starts"])
            origin = int(fossil["origin_time"])
            if is_candidate and origin >= cur_tail_origin:
                tail = fossil["id"]
        if tail == None:
            print("Failed to find representative organism in phylogeny.")
            exit(-1)
        # Build lineage from tail
        lineage = []
        current_id = tail
        while current_id != "NONE":
            lineage.append(current_id)
            current_id = phylogeny[current_id]["ancestor_list"].strip("[]").split(",")[0]
            if current_id == "NONE":
                phylogeny[lineage[-1]]["origin_time"] = -1

        print(lineage)

        # Collect summary statistics
        lineage_summary_stats = {"mrca":None,
                                 "length":0,
                                 "accum_muts":0,
                                 "accum_gene_move_muts":0,
                                 "accum_bit_flip_muts":0,
                                 "accum_bit_ins_muts":0,
                                 "accum_bit_del_muts":0
                                }
        lineage_summary_stats["length"] = len(lineage)
        for ancestor_id in lineage:
            # Collect mutation information
            gene_move_muts = int(phylogeny[ancestor_id]["gene_move_muts"])
            bit_flip_muts = int(phylogeny[ancestor_id]["bit_flip_muts"])
            bit_ins_muts = int(phylogeny[ancestor_id]["bit_ins_muts"])
            bit_del_muts = int(phylogeny[ancestor_id]["bit_del_muts"])
            total_muts = gene_move_muts + bit_flip_muts + bit_ins_muts + bit_del_muts
            # Add mutation information to summary stats
            lineage_summary_stats["accum_muts"] += total_muts
            lineage_summary_stats["accum_gene_move_muts"] += gene_move_muts
            lineage_summary_stats["accum_bit_flip_muts"] += bit_flip_muts
            lineage_summary_stats["accum_bit_ins_muts"] += bit_ins_muts
            lineage_summary_stats["accum_bit_del_muts"] += bit_del_muts
            # Add cumulative mutations to phylogeny
            phylogeny[ancestor_id]["accum_muts"] =  lineage_summary_stats["accum_muts"]
            phylogeny[ancestor_id]["accum_gene_move_muts"] =  lineage_summary_stats["accum_gene_move_muts"]
            phylogeny[ancestor_id]["accum_bit_flip_muts"] = lineage_summary_stats["accum_bit_flip_muts"]
            phylogeny[ancestor_id]["accum_bit_ins_muts"] = lineage_summary_stats["accum_bit_ins_muts"]
            phylogeny[ancestor_id]["accum_bit_del_muts"] = lineage_summary_stats["accum_bit_del_muts"]

        lineage.reverse()
        final_gen = update
        expanded_lineage = []
        for i in range(len(lineage)):
            ancestor_id = lineage[i]
            # origin_time = int(phylogeny[ancestor_id]["origin_time"]) + 1 # Shift everything by 1 generation to account for root starting at -1
            next_gen = int(phylogeny[lineage[i+1]]["origin_time"]) + 1 if i+1 < len(lineage) else final_gen
            while len(expanded_lineage) < next_gen:
                expanded_lineage.append(ancestor_id)

        # print("===========")
        # print(expanded_lineage)
        # print(len(expanded_lineage))
        # print("===========")


        # Fix gene starts for output format
        if "gene_starts" in rep_org:
            rep_org["gene_starts"] = "\"" + rep_org["gene_starts"] + "\""
        for fossil in phylogeny:
            phylogeny[fossil]["gene_starts"] = "\"" + phylogeny[fossil]["gene_starts"] + "\""

        # Collect lineage summary fields
        config_fields = [field for field in run_settings if (field not in config_exclude)]
        summary_fields = ["length","accum_muts","accum_gene_move_muts","accum_bit_flip_muts","accum_bit_ins_muts","accum_bit_del_muts"]
        lineage_summary_head_set.add(",".join(summary_fields + config_fields))
        summary_info =  [str(lineage_summary_stats[field]) for field in summary_fields] + [run_settings[field].strip() for field in config_fields]
        lineage_summary_info.append(",".join(summary_info))

        # Collect sequence fields
        # TODO!

        # Build a joint header.
        # # - gene stats fields
        # gene_stats_fields = [field for field in gene_stats_header if field not in gene_stats_exclude]
        # field_set = set(gene_stats_fields)
        # # - rep org fields
        # rep_org_fields = [field for field in rep_org_header if (field not in rep_org_exclude) and (field not in field_set)]
        # field_set.union(set(rep_org_fields))
        # # - run configuration fields
        # config_fields = [field for field in run_settings if (field not in config_exclude) and (field not in field_set)]
        # fields = gene_stats_fields + rep_org_fields + config_fields
        # # Combine all fields into single header (double check that this matches all previous computed headers)
        # combined_header_set.add(",".join(fields))
        # if len(combined_header_set) > 1:
        #     print("Header mismatch!")
        #     exit(-1)

    # Output aggregate information.
    # out_content = list(combined_header_set)[0] + "\n"
    # out_content += "\n".join([",".join(map(str, line)) for line in aggregate_info])
    # out_path = os.path.join(dump_dir, "agg_data.csv")
    # with open(out_path, "w") as fp:
    #     fp.write(out_content)
    # print(f"Done! Output written to {out_path}")

    summary_out_content = list(lineage_summary_head_set)[0] + "\n"
    summary_out_content += "\n".join(lineage_summary_info)
    summary_output_path = os.path.join(dump_dir, "lineages_summary.csv")
    with open(summary_output_path, "w") as fp:
        fp.write(summary_out_content)
    print(f"Done! Output written to {summary_output_path}")

if __name__ == "__main__":
    main()