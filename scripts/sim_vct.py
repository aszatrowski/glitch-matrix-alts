# Redirect stdout/stderr to log file
import sys
log_file = snakemake.log[0]
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

import xftsim as xft
from xftsim.reproduce import RecombinationMap
from xftsim.sim import Simulation # import as own object

import numpy as np
import random

from minimal_simulation_utils import make_covariance_matrix, build_minimal_founders

random.seed(1)
# import simulation parameters from snakemake environment object; make sure all are float or int
n_indivs=int(snakemake.params['n_indivs'])
n_causal = int(snakemake.params['n_causal'])
m_variants = int(snakemake.params['m_variants'])

h2=float(snakemake.wildcards['h2'])
b2=float(snakemake.wildcards['b2'])

print(f"n_indivs={n_indivs}, n_causal={n_causal}, m_variants={m_variants}, h2={h2}, b2={b2}")

founders, recomb = build_minimal_founders(
    n_indivs=n_indivs,
    m_variants=m_variants,
    min_af=0.1,
    max_af=0.5,
    chrom=1
)

def additive_effects_freq(founder_haplotypes, n_causal = None, w = 1, h2 = 0):
    """
    Minimal additive effects builder based on allele frequencies.
    """
    # Extract per-variant allele frequencies and minor allele frequencies.
    m =  founder_haplotypes.xft.m
    allele_freqs = founder_haplotypes.af
    m_allele_freqs = np.minimum(allele_freqs, 1-allele_freqs)[::2]
    # choose n_causal sites to have an effect and the rest to be zero
    beta = np.zeros((m,1))
    print(f"Number of causal sites: {n_causal}.")

    if n_causal is None or n_causal >= m:
        raise ValueError("n_causal must be set and less than m")

    # Apply thinning if n_causal is much smaller than m - just to speed up.
    if n_causal < m // 10:
        thinning_fraction = 0.2 
        thin_m = int(thinning_fraction * m)
        subset = np.random.choice(np.arange(m), thin_m, replace=False)
        subset_maf = m_allele_freqs[subset]  
        weights = (1 - subset_maf) ** w
        p = weights / weights.sum()
        causal_sites = np.random.choice(subset, n_causal, replace=False, p=p)
    else:
        weights = (1 - m_allele_freqs) ** w # FIXED: used to be subset_maf instead of m_allele_freqs, which is only assigned in the if: block
        p = weights / weights.sum()
        causal_sites = np.random.choice(np.arange(m), n_causal, replace=False, p=p)

    # Assign effect sizes to the chosen causal sites.
    beta[causal_sites] = np.random.randn(n_causal,1) * np.sqrt(h2)

    print(f"Mean allele frequency: {np.mean(m_allele_freqs)}.")
    print(f"Mean allele frequency at causal sites: {np.mean(m_allele_freqs[causal_sites])}.")

    # Package effects in xftsim structures.
    vindex = founder_haplotypes.xft.get_variant_indexer()
    cindex = xft.index.ComponentIndex.from_product('Y', 'additiveGenetic')
    effects = xft.effect.AdditiveEffects(beta=beta,
                                         variant_indexer=vindex,
                                         component_indexer=cindex)
    
    print(f"Done deciding effect sizes.")
    return([effects,causal_sites])

def vct_architecture(founder_haplotypes, h2, b2, parental_coeff = 1/2, n=1000, n_causal = None, w=0):
    """
    Function to set up VCT architecture. source: Marida, entirely.
    """
    # sets a dominance coefficient. b2m = coeff * b2m, b2m + b2p = b2
    b2p = parental_coeff * b2
    b2m = b2 - b2p
    if parental_coeff!=1/2:
        print(f"Paternal variance is {b2p} and maternal variance is {b2m}.")

    # Make genetic component.
    effects,causal_sites = additive_effects_freq(founder_haplotypes, n_causal = n_causal, h2 = h2, w=w)
    print("Now making architecture...")
    a_comp = xft.arch.AdditiveGeneticComponent(beta=effects,  component_name="additiveGenetic")
    # making the vertical component
    vert_input = xft.index.ComponentIndex.from_product(['HeritableFactor'], ['phenotype'], [0,1])
    vert_input.comp_type ='intermediate'
    vert_output = xft.index.ComponentIndex.from_product(['HeritableFactor'], ['vertical'], [-1])
    # one part will be from parents
    vt_comp = xft.arch.LinearVerticalComponent(input_cindex=vert_input,
                                              output_cindex=vert_output,
                                              founder_variances=[b2,b2],
                                              coefficient_matrix=np.array([[1-parental_coeff, parental_coeff]]),
                                              normalize = False)
    # also has variance
    vt_noise = xft.arch.AdditiveNoiseComponent(variances=[b2/2 * (1 + 1/n)], # extra factor to account for shrinking variance due to relatedness in sample
                                               phenotype_name=['HeritableFactor'],
                                               component_name='additiveNoise')
    vind = xft.index.ComponentIndex.from_product(['HeritableFactor'],
                                                 ['additiveNoise','vertical'])
    # makes an intermediate phenotype for the vertical component
    vtrans = xft.arch.SumAllTransformation(input_cindex=vind,output_component_name='phenotype',output_comp_type = 'intermediate')
    # additive (non heritable) noise
    a_noise = xft.arch.AdditiveNoiseComponent(variances=[1-h2-b2], phenotype_name=['Y'],component_name='additiveNoise')
    # constructing the final phenotype
    input_ind = xft.index.ComponentIndex(['HeritableFactor'], ['phenotype'],
                                         comp_type='intermediate')
    output_ind = xft.index.ComponentIndex(['Y'], ['heritableNoise'])
    # feed heritable factor into the phenotype stack
    coefficient_matrix = np.array([[1]])
    causal_comp = xft.arch.LinearTransformationComponent(input_ind, output_ind,
                                                   coefficient_matrix,normalize=False)
    input_cindex = xft.index.ComponentIndex(phenotype_name = ['Y', 'Y','Y'],
                                            component_name = ['heritableNoise',
                                                              'additiveNoise',
                                                              'additiveGenetic'])
    strans = xft.arch.SumAllTransformation(input_cindex,output_component_name='phenotype',
                                           output_comp_type = 'outcome')
    arch = xft.arch.Architecture([vt_noise,
                                    vt_comp,
                                    vtrans,
                                    a_noise,
                                    a_comp,
                                    causal_comp,
                                    strans])
    print("Done making architecture.")
    return([arch,causal_sites])

if h2 + b2 > 1.0:
    raise ValueError("Require h2 + b2 <= 1.0.")

# Pick a small causal set for this minimal example.

# if n_causal is None:
#     n_causal = max(5, m_variants // 20)
# if n_causal >= max(1, m_variants // 10):
#     n_causal = max(5, m_variants // 20)
#     print(f"[note] lowering n_causal to {n_causal} for this minimal demo.")


arch, _ = vct_architecture(
    founders,
    h2=h2,
    b2=b2,
    parental_coeff=0.5,
    n=n_indivs,
    n_causal=n_causal,
    w=0,
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
sim.run(5)
print('Complete.')
print(sim)
    
print('Building covariance matrix...')
cov = make_covariance_matrix(sim, maf = 0.01, include_pedigree=False)
print('Complete.')
print('Writing covariance matrix...')
cov.to_csv(snakemake.output['covariances_csv'], index=False)
print('Complete.')
