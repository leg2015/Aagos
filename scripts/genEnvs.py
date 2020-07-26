'''
Generate a specified number of random GRADIENT and NK environments.
'''

import argparse, random, os, errno

alphabet=[0,1]

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
    parser = argparse.ArgumentParser(description="Generate all genetic architectures")
    parser.add_argument("--num_genes", type=int, help="# of genes")
    parser.add_argument("--gene_size", type=int, help="# of bits in each genes")
    parser.add_argument("--envs", type=int, help="# of environments")
    parser.add_argument("--dump", type=str, help="Where to dump these environments?", default=".")

    args = parser.parse_args()
    num_genes = args.num_genes
    gene_bits = args.gene_size
    num_envs = args.envs
    dump = args.dump
    mkdir_p(dump)

    # Generate 'num_envs' number of unique NK environments
    gene_table_size = 2**gene_bits
    nk_environments = set()
    while len(nk_environments) < num_envs:
        nk_environments.add("[ " + " ".join([ "[ " + " ".join([str(random.random()) for _ in range(gene_table_size)]) + " ]" for gene in range(num_genes)]) + " ]")
    nk_environments = list(nk_environments)

    # Output nk environments
    nk_dump = os.path.join(dump, "nk")
    mkdir_p(nk_dump)
    for i in range(num_envs):
        fname = f"nk_env_{i}.env"
        with open(os.path.join(nk_dump, fname), "w") as fp:
            fp.write(nk_environments[i])
    nk_environments = None  # drop all that nk environment memory

    # Generate 'num_envs' number of unique GRADIENT environments
    gradient_environments = set()
    while len(gradient_environments) < num_envs:
        gradient_environments.add( "[ " + " ".join([ "".join([str(random.choice(alphabet)) for bit in range(gene_bits)]) for gene in range(num_genes) ]) + " ]" )
    gradient_environments = list(gradient_environments)

    # Output gradient environments
    gradient_dump = os.path.join(dump, "gradient")
    mkdir_p(gradient_dump)
    for i in range(num_envs):
        fname = f"gradient_env_{i}.env"
        with open(os.path.join(gradient_dump, fname), "w") as fp:
            fp.write(gradient_environments[i])

if __name__ == "__main__":
    main()
