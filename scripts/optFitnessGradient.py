

'''
Architecture given by
- genome length
- gene length
- gene count
- gene start positions

Functions:
- Given environment + architecture, compute optimal fitness
- Given architecture, compute expected fitness optimum

'''

import random, statistics, math


class GeneticArchitecture:
    def __init__(self,
                 genome_length,
                 gene_count,
                 gene_length,
                 gene_starts):
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
        print(f"Expected optimal fitness: {self.ComputeExpectedOptimalFitness()}")

    def ComputeExpectedOptimalFitness(self):
        return sum(self._computeExpectedSiteContribution(self.coding_site_occupant_cnt[site]) for site in self.coding_site_occupant_cnt)

    def _computeExpectedSiteContribution(self, num_occupants):

        if num_occupants==0: return 0.0   # If there are no occupants, it cannot contribute anything
        if num_occupants==1: return 1.0   # If there's only one occupant, it will contribute the full 1.0

        # There is more than one occupant at this site...
        # first, compute the distribution of possible gene agreements
        possible_configurations = 2**num_occupants     # this doesn't work for bigger alphabet sizes

        worst_agreement = math.ceil(num_occupants/2.0) # worst agreement is half (rounded up) of bits agreeing

        def num_agreements(shared):
            return math.factorial(num_occupants) / float(math.factorial(shared) * math.factorial(num_occupants - shared))

        # possible_agreements = {agreement:0 for agreement in range(worst_agreement, num_occupants+1)}
        possible_agreements=None
        if num_occupants & 1:
            # odd-length bit string
            possible_agreements={num_shared: 2*num_agreements(num_shared) for num_shared in range(worst_agreement, num_occupants+1)}
        else:
            # even-length bit string
            possible_agreements={num_shared: (1 if num_shared == worst_agreement else 2) * num_agreements(num_shared) for num_shared in range(worst_agreement, num_occupants+1)}

        return sum([possible_agreements[num_shared]*num_shared for num_shared in possible_agreements]) / possible_configurations



    def ComputeOptimalFitness(self, gene_targets):
        # gene_targets = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        possible_position_values = list({val for target in gene_targets for val in target})
        possible_position_values.sort()
        # print(f"Possible target values: {possible_position_values}")

        # Optimal genome (coding region): for each coding site, set position equal to majority target value
        target_positions_satisfied = 0 # How many positions (across all targets) can we satisfy?

        optimal_site_map={pos:"X" for pos in range(self.genome_length)}

        for coding_pos in self.coding_site_occupancy:
            # print(f"--computing coding pos {coding_pos}--")
            position_votes = {val:0 for val in possible_position_values}
            for gene in self.coding_site_occupancy[coding_pos]:
                gene_id = gene["gene_id"]
                gene_index = gene["gene_index"]

                # What does this gene want this position to be?
                gene_vote = gene_targets[gene_id][gene_index]
                position_votes[gene_vote] += 1
            # print(f"Coding position: {coding_pos}; Position votes: {position_votes}")
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

        # print(f"Optimal site map: {optimal_site_map}")
        # print(f"Optimal fitness: {target_positions_satisfied}")
        return {"optimal_sites": optimal_site_map, "optimal_fitness": target_positions_satisfied}

def main():

    genome_length = 16
    gene_count=16
    gene_length=4
    alphabet=[0,1]
    # starts = [ [0,0,0,0], [0,0,0,1], [0,0,0,2], [0,0,2,2], [0,4,8,12], [0,1,5,7] ]
    # starts = [ [0,0], [0,2], [0,4], [0,5], [0,6], [0,7], [0,8] ]
    starts = [ [random.randint(0,genome_length) for g in range(gene_count)] for _ in range(50)]

    num_envs = 1000

    architectures = [GeneticArchitecture(genome_length,gene_count,gene_length,start) for start in starts]

    # random gene targets
    gene_targets = [[[random.choice(alphabet) for _ in range(gene_length)] for _ in range(gene_count)] for i in range(num_envs)]
    for architecture in architectures:
        print(f"===== architecture: {architecture.gene_starts} =====")
        # Compute expected fitness for achitecture
        expected_optimal_fitness = architecture.ComputeExpectedOptimalFitness()
        # architecture.Print()
        # Compute optimal on set of randomly generated gene targets
        fitness_distribution = [architecture.ComputeOptimalFitness(gene_targets[i])["optimal_fitness"] for i in range(num_envs)]
        mean = statistics.mean(fitness_distribution)
        # median = statistics.median(fitness_distribution)
        # mode = statistics.mode(fitness_distribution)
        conf = 1.96 * ( statistics.stdev(fitness_distribution) / math.sqrt(num_envs) )
        print(f"Expected={expected_optimal_fitness}")
        print(f"mean    ={mean}  [{mean-conf}:{mean+conf}]")


if __name__ == "__main__":
    main()