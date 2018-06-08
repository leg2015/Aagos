#include <iostream>

#include "base/vector.h"
#include "config/ArgManager.h"
#include "config/command_line.h"

#include "../AagosOrg.h"
#include "../AagosWorld.h"

int main(int argc, char* argv[])
{
  AagosConfig config;
  // Deal with loading config values via native interface (config file and command line args)
  config.Read("Aagos.cfg");
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "Aagos.cfg", "Aagos-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  AagosWorld world(config);

  emp::Random & random = world.GetRandom();

  // Build a random initial population
  for (uint32_t i = 0; i < config.POP_SIZE(); i++) {
    AagosOrg next_org(config.NUM_BITS(), config.NUM_GENES(), config.GENE_SIZE()); // build org
    next_org.Randomize(random); // randomize org
    world.Inject(next_org);     // inject org
  }

  // runs ech generation
  for (size_t gen = 0; gen < config.MAX_GENS(); gen++) {
    // Do mutations on the population.
    world.DoMutations(config.ELITE_COUNT());

    // TODO: looks like histogram is still acting up, need to get it working tomorrow so can start runs]
    // const emp::vector<size_t> & hist = world.GetOrg(1).GetHistogram().GetHistCounts();
    // std::cout << "bin count in hist bin 0: " << world.GetOrg(1).GetHistogram().GetHistCount(0) << std::endl;
    // std::cout << "width of hist bin 0: " << world.GetOrg(1).GetHistogram().GetHistWidth(1) << std::endl;
    // std::cout << "min possible hist val: " << world.GetOrg(1).GetHistogram().GetHistMin() << std::endl;
    // std::cout << "max possible hist val: " << world.GetOrg(1).GetHistogram().GetHistMax() << std::endl;
    // Keep the best individual.

    if (config.ELITE_COUNT()) emp::EliteSelect(world, config.ELITE_COUNT(), 1);

    // Run a tournament for the rest...
    emp::TournamentSelect(world, config.TOURNAMENT_SIZE(), config.POP_SIZE()-config.ELITE_COUNT());

    // Update world
    world.Update();

    // If it's a generation to print to console, do so
    if (gen % config.PRINT_INTERVAL() == 0) {
      std:: cout << "-----------gen" << gen << "----------------" << std::endl;
       std::cout << "gene neighbors: "<< emp::to_string(world.GetOrg(0).GetGeneNeighbors()) << std::endl;
      for(size_t i = 0; i < config.POP_SIZE(); i++) {
      std::cout << gen
                << " : fitness=" << world.CalcFitnessID(i)
                << " size=" << world[i].GetNumBits() << std::endl;
                // << "histogram count=[";
                // for(int i = 0; i < hist.size(); i++) {
                //  std::cout << hist[i] << ", ";
                // }
                // std::cout  << "]" << std::endl;
      world[i].Print();
      }
    }
  }
}
