import xftsim as xft
from xftsim.reproduce import RecombinationMap
from xftsim.sim import Simulation # import as own object

import numpy as np
import random

from minimal_simulation_utils import make_covariance_matrix, build_minimal_founders

random.seed(1)
founders, recomb = build_minimal_founders(
    n_indivs=int(snakemake.params['n_indivs']),
    m_variants=int(snakemake.params['m_variants']),
    min_af=0.1,
    max_af=0.5,
    chrom=1
)
# build minimal additive genetic + additive noise architecture
arch = xft.arch.GCTA_Architecture(
    h2=float(snakemake.wildcards['h2']),
    phenotype_name='Y',
    haplotypes=founders
)

mating_regime = xft.mate.RandomMatingRegime(
    mates_per_female=1,
    offspring_per_pair=2,
    female_offspring_per_pair="balanced",
    exhaustive=True,
)

post_processors = [xft.proc.LimitMemory(n_haplotype_generations=1)]
sim = xft.sim.Simulation(
    founder_haplotypes=founders,
    architecture=arch,
    recombination_map=recomb,
    mating_regime=mating_regime,
    statistics=[],
    post_processors=post_processors
)

print('Running generations...')
sim.run(10)
print('Complete.')
    
print('Building covariance matrix...')
cov = make_covariance_matrix(sim, maf = 0.01, include_pedigree=False)
print('Complete.')
print('Writing covariance matrix...')
cov.to_csv(snakemake.output['covariances_csv'], index=False)
print('Complete.')
