
Test robustness of results at higher/lower gene move and insertion/deletion mutation rates.

- Bit flip rates
  - `-BIT_FLIP_PROB 0.003 -PHASE_2_BIT_FLIP_PROB 0.003`
  - `-BIT_FLIP_PROB 0.03 -PHASE_2_BIT_FLIP_PROB 0.03`
- Gene move rates
  - `-GENE_MOVE_PROB 0.0003`
  - `-GENE_MOVE_PROB 0.003`
  - `-GENE_MOVE_PROB 0.03`
- Insertion/deletion rates
  - `-BIT_INS_PROB 0.0001 -BIT_DEL_PROB 0.0001"`
  - `-BIT_INS_PROB 0.001 -BIT_DEL_PROB 0.001"`
  - `-BIT_INS_PROB 0.01 -BIT_DEL_PROB 0.01`