import xftsim as xft
from xftsim.reproduce import RecombinationMap
from xftsim.sim import Simulation # import as own object

import numpy as np
import matplotlib.pyplot as plt

from minimal_simulation_utils import make_covariance_matrix

def build_minimal_founders(n_indivs: int, m_variants: int, min_af: float = 0.05, max_af: float = 0.5, chrom: int = 1):
    """
    Create founder haplotypes with uniform positions and a simple recomb map.
    """
    afs = np.random.uniform(min_af, max_af, size = m_variants)
    founders = xft.founders.founder_haplotypes_from_AFs(n=n_indivs, afs=afs)

    # 2 * m_variants
    pos_len = founders["pos_bp"].shape[0]
    # assign genomic positions
    pos_bp = (np.arange(1, pos_len + 1) * 1000).astype(int)
    # 1 cM (1/100 recombinations) for every 1M variants
    pos_cM = pos_bp / 1_000_000

    # assign positions and cM rates to founders
    founders["pos_bp"].values = pos_bp
    founders["pos_cM"].values = pos_cM
    # repeat across chromosomes
    founders["chrom"].values = np.repeat(chrom, pos_len)

    # uniform recombination map consistent with these parameters.
    # recomb is not immediately necessary, but is required for the actual xft.sim.Simulation
    recomb = RecombinationMap.variable_map_from_haplotypes_with_cM(founders)
    return founders, recomb

founders, recomb = build_minimal_founders(
    n_indivs=100,
    m_variants=1_000,
    min_af=0.1,
    max_af=0.5,
    chrom=1
)
# build minimal additive genetic + additive noise architecture
arch = xft.arch.GCTA_Architecture(
    h2=0.5,
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

sim.run(4)

cov = make_covariance_matrix(sim, maf = 0.01, include_pedigree=False)
cov.to_csv("data/gcta_covariances.csv", index=False)
