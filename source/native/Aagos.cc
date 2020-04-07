#include <iostream>

#include "base/vector.h"
#include "config/ArgManager.h"
#include "config/command_line.h"

#include "../AagosConfig.h"
#include "../AagosOrg.h"
#include "../AagosWorld.h"

int main(int argc, char* argv[])
{
  std::string config_fname = "Aagos.cfg";
  AagosConfig config;
  // Deal with loading config values via native interface (config file and command line args)
  config.Read(config_fname);
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "Aagos.cfg", "Aagos-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  // Write to screen how the experiment is configured
  std::cout << "==============================" << std::endl;
  std::cout << "|    How am I configured?    |" << std::endl;
  std::cout << "==============================" << std::endl;
  config.Write(std::cout);
  std::cout << "==============================\n" << std::endl;

  emp::Random random(config.SEED());
  AagosWorld world(random, config);

  // Build a random initial population
  for (uint32_t i = 0; i < config.POP_SIZE(); i++) {
    AagosOrg next_org(config.NUM_BITS(), config.NUM_GENES(), config.GENE_SIZE()); // build org
    next_org.Randomize(random); // randomize org
    world.Inject(next_org);     // inject org
  }

  // Configure world to automatically mutate offspring on placement.
  world.SetAutoMutate(config.ELITE_COUNT());

  // runs each generation
  for (size_t gen = 0; gen <= config.MAX_GENS(); gen++) {
    // TODO - pre-calculate fitness for each organism

    // ---- do selection ----
    // Keep the best individual.
    if (config.ELITE_COUNT()) emp::EliteSelect(world, config.ELITE_COUNT(), 1);
    // Run a tournament for the rest...
    emp::TournamentSelect(world, config.TOURNAMENT_SIZE(), config.POP_SIZE()-config.ELITE_COUNT());

    // Update world (handles mutations)
    world.Update();

    // If it's a generation to print to console, do so
    if (gen % config.PRINT_INTERVAL() == 0) {
      std::cout << gen
                << " : fitness=" << world.CalcFitnessID(0)
                << " size=" << world[0].GetNumBits()
                << std::endl;
      world[0].Print();
    }

  }
}
