# data_higher_divergence

This folder contains repeat-level summary tables for the divergence-gradient simulation set. 

## exNG input files

- `small_exNG_summary_0.05_0.1_0.5.tsv`: 39 small-species IDs at divergences `0.05`, `0.1`, and `0.5` (and legacy divergence-`1` rows, which are excluded).
- `small_exNG_summary_0.2_0.3_0.4.tsv`: the same 39 small-species IDs at divergences `0.2`, `0.3`, and `0.4`.
- `large_exNG_0.05_0.1_0.2_0.3_0.4_0.5.tsv`: the three additional large species, completing the 42-species analysis, for all six divergences.

small species and large species are listed in `data/species_list_small` and `data/species_list_large`. `data/species_list_all` lists all of them.

Species are identified by csv row indices in `data/summary_of_mutation_spectrum_and_omega_ratio.csv`.

Small and large species are analyzed together although their data are in different files.

A combination is one `species + omega + divergence` configuration with 100 repeats.

The divergence gradient is `0.05`, `0.1`, `0.2`, `0.3`, `0.4`, and `0.5`. 

Required analysis columns are `species`, `divergence`, `omega`, `selection`, `repeat`, `pNpS`, and `dNdS`.

