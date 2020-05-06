//  This file is part of Project Name
//  Copyright (C) Michigan State University, 2017.
//  Released under the MIT Software license; see doc/LICENSE

#include "web/web.h"
#include "../AagosConfig.h"
#include "../AagosWebD3.h"
// #include "Evolve/NK.h"

struct AagosWebWrapper {
  AagosConfig cfg;
  emp::Ptr<AagosWebInterface> interface;

  AagosWebWrapper() {
    // Configure aagos
    cfg.POP_SIZE(100);
    cfg.MAX_GENS(1000);
    cfg.SEED(2);
    cfg.TOURNAMENT_SIZE(8);
    cfg.GRADIENT_MODEL(true);
    cfg.CHANGE_MAGNITUDE(1);
    cfg.CHANGE_FREQUENCY(8);
    cfg.NUM_BITS(32);
    cfg.NUM_GENES(4);
    cfg.GENE_SIZE(8);
    cfg.MAX_SIZE(128);
    cfg.MIN_SIZE(4);
    cfg.GENE_MOVE_PROB(0.003);
    cfg.BIT_FLIP_PROB(0.003);
    cfg.BIT_INS_PROB(0.001);
    cfg.BIT_DEL_PROB(0.001);
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
