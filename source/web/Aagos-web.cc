//  This file is part of Project Name
//  Copyright (C) Michigan State University, 2017.
//  Released under the MIT Software license; see doc/LICENSE

#include "web/web.h"
#include "../AagosConfig.h"
#include "../AagosWeb.h"
// #include "Evolve/NK.h"

struct AagosWebWrapper {
  AagosConfig cfg;
  emp::Ptr<AagosWebInterface> interface;

  AagosWebWrapper() {
    // Configure aagos
    cfg.POP_SIZE(500);  // original: 100
    cfg.MAX_GENS(1500); // original: 1000
    cfg.SEED(42);
    cfg.TOURNAMENT_SIZE(8);
    cfg.GRADIENT_MODEL(true);
    cfg.CHANGE_MAGNITUDE(1);
    cfg.CHANGE_FREQUENCY(8);
    cfg.NUM_BITS(64);         // original: 32
    cfg.NUM_GENES(8);         // original: 4
    cfg.GENE_SIZE(8);
    cfg.MIN_SIZE(8);
    cfg.MAX_SIZE(512);        // original: 128
    cfg.GENE_MOVE_PROB(0.001); // original: 0.003
    cfg.BIT_FLIP_PROB(0.001);  // original: 0.003
    cfg.BIT_INS_PROB(0.001);
    cfg.BIT_DEL_PROB(0.001);
    cfg.PHASE_2_ACTIVE(false);
    cfg.PHASE_2_LOAD_ENV_FROM_FILE(false);
    cfg.PHASE_2_CHANGE_MAGNITUDE(0);
    cfg.PHASE_2_CHANGE_FREQUENCY(0);
    cfg.PHASE_2_MAX_GENS(1000);
    cfg.PHASE_2_TOURNAMENT_SIZE(cfg.TOURNAMENT_SIZE());
    cfg.PHASE_2_GENE_MOVE_PROB(0.0);
    cfg.PHASE_2_BIT_FLIP_PROB(cfg.BIT_FLIP_PROB());
    cfg.PHASE_2_BIT_INS_PROB(0.0);
    cfg.PHASE_2_BIT_DEL_PROB(0.0);
    interface = emp::NewPtr<AagosWebInterface>(cfg);
  }

  ~AagosWebWrapper() {
    interface.Delete();
  }
};

AagosWebWrapper aagos;

int main() {
  std::cout << "Hello?" << std::endl;
}
