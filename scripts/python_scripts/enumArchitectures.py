
import itertools, argparse, os, sys, errno, csv
import networkx as nx

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
Compact gene starting positions
'''
def compact_positions(gene_bits, num_genes, gene_starts):
    gene_starts = sorted(list(gene_starts))
    prev_end = 0
    for i in range(len(gene_starts)):
        max_pos = min(i * gene_bits, prev_end)
        gene_starts[i] = max_pos if gene_starts[i] > max_pos else gene_starts[i]
        prev_end = gene_starts[i] + gene_bits
    return gene_starts

def positions_to_graph(gene_starts, gene_bits, genome_length):
    gene_positions = {gene:{ (gene_starts[gene] + i) % genome_length for i in range(gene_bits) } for gene in range(len(gene_starts))}
    site_occupancy = {p:{gene_id for gene_id in range(len(gene_starts)) if p in gene_positions[gene_id] } for p in range(genome_length)}
    edges = {}
    for site in site_occupancy:
        occupants = list(site_occupancy[site])
        for oi in range(0, len(occupants)):
            for ok in range(oi+1, len(occupants)):
                edge = frozenset((occupants[oi], occupants[ok]))
                if edge not in edges: edges[edge] = 0
                edges[edge] += 1

    gene_graph = nx.Graph()
    gene_graph.add_nodes_from([i for i in range(len(gene_starts))]) # Add graph nodes
    for edge in edges:
        u, v = tuple(edge)
        w = edges[edge]
        gene_graph.add_edge(u, v, weight = w)

    # print(f"Gene position: {gene_positions}")
    # print(f"Site occupancy: {site_occupancy}")
    # print(f"Overlap edges: {edges}")
    # print(f"Gene graph nodes: {str(gene_graph.nodes)}")
    # print(f"Gene graph adj: {str(gene_graph.adj)}")

    return gene_graph

def is_isomorphic(g1, g2):
    if not nx.faster_could_be_isomorphic(g1, g2): return False
    if not nx.fast_could_be_isomorphic(g1, g2): return False
    if not nx.could_be_isomorphic(g1, g2): return False
    return nx.is_isomorphic(g1, g2, edge_match=nx.algorithms.isomorphism.numerical_edge_match('weight', 1))

def main():
    parser = argparse.ArgumentParser(description="Generate all genetic architectures")
    parser.add_argument("--num_genes", type=int, help="# of genes")
    parser.add_argument("--gene_size", type=int, help="# of bits in each genes")
    parser.add_argument("--genome_size", type=int, help="# of bits in genome")
    parser.add_argument("--output_as_genomes", action="store_true", help="Should we output results as genome?")
    parser.add_argument("--output_as_gradient_targets", action="store_true", help="Should we output results as gradient gene targets?")
    parser.add_argument("--single_file", action="store_true", help="Should we output results as a single file?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")

    args = parser.parse_args()
    num_genes = args.num_genes
    gene_bits = args.gene_size
    genome_length = args.genome_size

    output_as_genomes = args.output_as_genomes
    output_as_gradient_targets = args.output_as_gradient_targets
    single_file = args.single_file


    dump_dir = args.dump

    # gene_bits = 4
    # num_genes = 8
    # genome_length = 32
    # max_size = gene_bits * num_genes

    # Compute all starting positions for a non-circular genome.
    all_starts = { positions for positions in itertools.combinations_with_replacement(range(genome_length), num_genes) }
    all_starts = sorted(list(all_starts))

    # don't print all starts if there are a ton...
    if len(all_starts) < 100000:
        print(f"All gene starts: {all_starts}")
    print(f"Number of possible gene starts: {len(all_starts)}")

    # For each gene layout, compute the associated graph.
    architectures = {}
    counter = 0
    total = len(all_starts)
    for layout in all_starts:
        print(f"{counter}/{total}")
        # Is this layout isomorphic with any of the architectures saved so far?
        layout_graph = positions_to_graph(layout, gene_bits, genome_length)
        # Actually, we want to bail as soon as we can
        is_iso = False
        for other in architectures:
            is_iso = is_isomorphic(layout_graph, architectures[other])
            if is_iso: break
        if not is_iso:
            architectures[layout] = layout_graph
        counter += 1

    # diff_arches = [layout for layout in architectures]
    print("===== Different architectures ====")
    for arch in architectures:
        print(f"Gene starts: {arch}")
        print(f"Graph adj: {architectures[arch].edges}")
        print("-----")
    print(f"Number of different architectures: {len(architectures)}")

    # Output architectures
    arch_bits = ''.join(['0' for i in range(genome_length)])
    lines = []
    for arch in architectures:
        gene_starts_str = ",".join(map(str, arch))
        lines.append(gene_starts_str + "," + arch_bits)

    mkdir_p(dump_dir)
    if output_as_genomes and single_file:
        out_path = os.path.join(dump_dir, "architectures.csv")
        with open(out_path, "w") as fp:
            fp.write("\n".join(lines))
    elif output_as_genomes and not single_file:
        for arch_i in range(len(lines)):
            out_path = os.path.join(dump_dir, f"architecture-{arch_i}.csv")
            with open(out_path, "w") as fp:
                fp.write(lines[arch_i])





    # print("===G1===")
    # g1 = positions_to_graph((0, 0, 4), 4, 12)
    # print("===G2===")
    # g2 = positions_to_graph((0, 4, 4), 4, 12)


    # is_iso = nx.is_isomorphic(g1, g2, edge_match=edge_match_fun)
    # print("=======")
    # print(f"G1: {g1.adj}")
    # print(f"G2: {g2.adj}")
    # print(f"is iso? {is_iso}")




    # print(all_starts)
    # print(len(all_starts))

    # [positions for positions in  ]

    # len({ tuple([v - i[0] for v in i ])  for i in itertools.combinations_with_replacement(range(16), 4) })



if __name__ == "__main__":
    main()