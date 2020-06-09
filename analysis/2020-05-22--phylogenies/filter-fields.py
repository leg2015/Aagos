import argparse, os, sys, errno, csv


include_fields = ["generation",
                  "coding_sites",
                  "genome_length",
                  "BIT_FLIP_PROB",
                  "TOURNAMENT_SIZE",
                  "CHANGE_MAGNITUDE",
                  "CHANGE_FREQUENCY",
                  "GRADIENT_MODEL",
                  "SEED"]

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
    parser = argparse.ArgumentParser(description="Data aggregation script")
    parser.add_argument("--data", type=str, help="Path to data file")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")

    args = parser.parse_args()
    data_path = args.data
    dump_dir = args.dump

    # Make place to dump aggregated data.
    mkdir_p(dump_dir)

    # Make filtered file stub
    filtered_header = ",".join(include_fields) + "\n"
    filtered_out_path = os.path.join(dump_dir, "filtered_lineage.csv")
    with open(filtered_out_path, "w") as fp:
        fp.write(filtered_header)

    print("Filtering data...")
    header = []
    filtered_line_buffer = []
    with open(data_path, "r") as fp:
        header = fp.readline().strip().split(",")
        header_lu = {header[i].strip():i for i in range(0, len(header))}
        for cnt, line in enumerate(fp):
            line = line.strip().split(",")
            filtered_line_buffer.append(",".join([line[header_lu[field]] for field in include_fields]))
            # If buffer gets too big, flush to out file
            if len(filtered_line_buffer) >= 1000:
                with open(filtered_out_path, "a") as out_fp:
                    out_fp.write("\n".join(filtered_line_buffer) + "\n")
                filtered_line_buffer = []

    with open(filtered_out_path, "a") as out_fp:
        out_fp.write("\n".join(filtered_line_buffer) + "\n")
    filtered_line_buffer = []

    print("Done!")



if __name__ == "__main__":
    main()