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

  AagosWorld world(config);
  world.Run();
}
