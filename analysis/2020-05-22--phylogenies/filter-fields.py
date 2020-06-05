import argparse, os, sys, errno, csv


include_fields = ["generation",
                  "coding_sites",
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

    # Read data
    content = None

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




    # print("...finished filtering data")

    # header = content[0].split(",")
    # # header_lu = {header[i].strip():i for i in range(0, len(header))}
    # content = content[1:]

    # include_fields_set = set(include_fields)
    # content = [ {header[i]: l[i] for i in range(0, len(l)) if header[i] in include_fields_set} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True) ]


    # print("building output file...")

    # filtered_lines = []
    # for line in content:
    #     filtered_lines.append(",".join([line[field] for field in include_fields]))
    # content = None
    # filtered_content = filtered_header + "\n"
    # filtered_content += "\n".join(filtered_lines)
    # print("writing output file...")
    # with open(os.path.join(dump_dir, "filtered_lineage.csv"), "w") as fp:
    #     fp.write(filtered_content)
    print("Done!")



if __name__ == "__main__":
    main()