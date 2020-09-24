"""
ASSUMES 2-GENE SYSTEM!
"""

import random, argparse

alphabet = [0, 1]
mutation_rates = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0]
overlap_rates = [0, 0.5, 1.0]

class GeneticArchitecture:
    def __init__(
        self,
        genome_length,
        gene_count,
        gene_length,
        gene_starts
    ):
        self.genome_length = genome_length
        self.gene_count = gene_count
        self.gene_length = gene_length
        self.gene_starts = gene_starts

        # What are each gene's genome positions
        self.gene_positions = {gene:[ (self.gene_starts[gene] + i) % self.genome_length for i in range(self.gene_length) ] for gene in range(len(self.gene_starts))}

        # Which sites in the genome are coding?
        self.coding_sites = {pos for gene in self.gene_positions for pos in self.gene_positions[gene]}

        # For each site in the genome, which genes occupy that site?
        self.site_occupancy = {pos:{gene for gene in self.gene_positions if pos in self.gene_positions[gene]} for pos in range(self.genome_length) }

        # For each coding site, which gene positions (for each occupying gene) map to that site?
        self.coding_site_occupancy = {pos:[] for pos in self.coding_sites}
        for pos in self.coding_site_occupancy:
            for gene_id in self.site_occupancy[pos]:
                gene_index = self.gene_positions[gene_id].index(pos) # By definition position should exist in this list.
                self.coding_site_occupancy[pos].append({"gene_id": gene_id, "gene_index": gene_index})

        # For each coding site, how many gene occupy that site?
        self.coding_site_occupant_cnt={pos:len(self.site_occupancy[pos]) for pos in self.coding_sites}

    def Print(self):
        print(f"Genome Length: {self.genome_length}; Gene Length: {self.gene_length}; Gene Count: {self.gene_count}")
        print(f"Gene start positions: {self.gene_starts}")
        print(f"Gene positions: {self.gene_positions}")
        print(f"Coding sites: {self.coding_sites}")
        print(f"Site occupancy: {self.site_occupancy}")
        print(f"Coding site occupancy map: {self.coding_site_occupancy}")
        print(f"Coding site occupant count: {self.coding_site_occupant_cnt}")

    def ComputeOptimal(self, gene_targets):
        # gene_targets = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        possible_position_values = list({val for target in gene_targets for val in target})
        possible_position_values.sort()

        # Optimal genome (coding region): for each coding site, set position equal to majority target value
        target_positions_satisfied = 0 # How many positions (across all targets) can we satisfy?

        optimal_site_map={pos:"X" for pos in range(self.genome_length)}

        for coding_pos in self.coding_site_occupancy:
            position_votes = {val:0 for val in possible_position_values}
            for gene in self.coding_site_occupancy[coding_pos]:
                gene_id = gene["gene_id"]
                gene_index = gene["gene_index"]
                # What does this gene want this position to be?
                gene_vote = gene_targets[gene_id][gene_index]
                position_votes[gene_vote] += 1
            # Which value should this coding site take on?
            best_value = None
            for val in possible_position_values:
                if best_value == None:
                    best_value = val
                elif position_votes[val] > position_votes[best_value]:
                    best_value = val
            optimal_site_map[coding_pos]=best_value
            # How many targets does this value satisfy (i.e., how many votes did this value get)?
            target_positions_satisfied += position_votes[best_value]
        return {"optimal_sites": optimal_site_map, "optimal_fitness": target_positions_satisfied}

    def ComputeFitness(self, genome, gene_targets):
        # For each gene target, check if site matches.
        fitness = 0
        for gene_i in range(len(gene_targets)):
            for ti in range(len(gene_targets[gene_i])):
                fitness += int(genome[self.gene_positions[gene_i][ti]] == gene_targets[gene_i][ti])
        return fitness

def RandomizeTargets(targets):
    """Randomizes list of gene targets.

    - targets: list of gene target (list of values in the alphabet)
    """
    for target in targets:
        for i in range(len(target)): target[i] = random.choice(alphabet)



def main():
    parser = argparse.ArgumentParser(description="Sample offspring fitness for a range of two gene architectures at a range of mutation rates.")
    parser.add_argument("--samples", type=int, help="How many independent samples to take per condition?")
    # parser.add_argument("--num_genes", type=int, help="# of genes")
    parser.add_argument("--gene_size", type=int, help="# of bits in each genes")
    # parser.add_argument("--genome_size", type=int, help="# of bits in genome")
    # parser.add_argument("--output_as_genomes", action="store_true", help="Should we output results as genome?")
    # parser.add_argument("--output_as_gradient_targets", action="store_true", help="Should we output results as gradient gene targets?")
    # parser.add_argument("--single_file", action="store_true", help="Should we output results as a single file?")
    # parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")

    args = parser.parse_args()
    samples = args.samples
    gene_size = args.gene_size

    # Compute the architectures
    overlap_amounts = {int(o * gene_size) for o in overlap_rates}
    architectures = {amount: [0, gene_size - amount] for amount in overlap_amounts}
    targets = [[alphabet[0] for _ in range(gene_size)] for ti in range(2)]

    fields = ["sample_id", "mutation_rate", "overlap", "parent_fitness", "offspring_fitness", "fitness_difference", "num_mutations"]
    data = []

    # For each architecture and mutation rate, draw N independent samples (environment + mutation).
    # - Assume optimal genotype given environment
    for overlap in architectures:
        print(f"----- Overlap: {overlap} ----- ")
        genome_size = (2*gene_size) - overlap
        architecture = GeneticArchitecture(
            genome_length = genome_size,
            gene_count = 2,
            gene_length = gene_size,
            gene_starts = architectures[overlap]
        )

        architectures[overlap]
        for mutation_rate in mutation_rates:
            for sample_i in range(samples):
                # architecture.Print()
                # (1) Randomize the environment.
                RandomizeTargets(targets)
                # (2) Compute optimal genome given environment and architecture.
                parent = architecture.ComputeOptimal(targets)
                # (3) Apply mutations to parent to generate offspring.
                offspring = [ int(not parent["optimal_sites"][site]) if random.random() < mutation_rate else parent["optimal_sites"][site] for site in range(genome_size)]
                # (4) Compute offspring fitness.
                offspring_fitness = architecture.ComputeFitness(offspring, targets)
                # (5) Record data.
                data.append({
                    "sample_id": sample_i,
                    "mutation_rate": mutation_rate,
                    "overlap": overlap,
                    "parent_fitness": parent["optimal_fitness"],
                    "offspring_fitness": offspring_fitness,
                    "fitness_difference": parent["optimal_fitness"] - offspring_fitness,
                    "num_mutations": sum(int(offspring[i] != parent["optimal_sites"][i]) for i in range(len(offspring)))
                })
                # print(f"Targets: {targets}")
                # print(f"Parent: {parent}")
                # print(f"Offspring: {offspring}")
                # print(f"Offspring fitness: {offspring_fitness}")
    # Write output file.
    content = ",".join(fields) + "\n"
    content += "\n".join([ ",".join([str(line[field]) for field in fields]) for line in data ])
    with open("mutant_samples.csv", "w") as fp:
        fp.write(content)


if __name__ == "__main__":
    main()