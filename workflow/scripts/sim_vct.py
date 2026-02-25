# Redirect stdout/stderr to log file
import sys
log_file = snakemake.log[0]
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

# suppress warnings about the impending deprecation of copy-in-place. not ideal, but actively encouraged on the xftsim README until fixed. see: https://github.com/border-lab/xftsim/blame/7aff3ea27932a032e96b923f310d247ff535028f/README.md. Note that his has been in place for ~2 years now.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import datetime
import xftsim as xft
import xarray as xr # for haplotype arrays

import numpy as np
import pandas as pd
import random

from minimal_simulation_utils import build_minimal_founders

arch = snakemake.wildcards['arch']
h2=float(snakemake.wildcards['h2'])
b2=float(snakemake.wildcards['b2'])
parental_coef = float(snakemake.wildcards['parental_coef'])
rep = snakemake.wildcards['rep']
random.seed(rep) # ensure each replicate is different yet individually reproducible

# import simulation parameters from snakemake environment object; make sure all are float or int
n_indivs=int(snakemake.params['n_indivs'])
n_causal = int(snakemake.params['n_causal'])
m_variants = int(snakemake.params['m_variants'])
chrom_count = int(snakemake.params['chrom_count'])

print(f"arch={arch}, n_indivs={n_indivs}, n_causal={n_causal}, m_variants={m_variants}, h2={h2}, b2={b2}, parental_coef={parental_coef}, rep/seed={rep}")

founders, recomb = build_minimal_founders(
    n_indivs=n_indivs,
    m_variants=m_variants,
    min_af=0.1,
    max_af=0.5,
    chrom=2
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

class my_RandomMatingRegime(xft.mate.MatingRegime):
    """
    SOURCE: Marida, glitch-matrix/workflow/scripts/simulation_utils.py
    A mating regime that randomly pairs individuals and produces offspring with balanced numbers of males and females.

    Parameters
    ----------
    offspring_per_pair : xft.utils.VariableCount, optional
        Number of offspring produced per mating pair, by default xft.utils.ConstantCount(1)
    mates_per_female : xft.utils.VariableCount, optional
        Number of males that mate with each female, by default xft.utils.ConstantCount(1)
    female_offspring_per_pair : Union[str, xft.utils.VariableCount], optional
        The number of female offspring per mating pair. If "balanced", the number is balanced with
        the number of male offspring. By default, "balanced".
    sex_aware : bool, optional
        If True, randomly paired individuals are selected so that there is an equal number of males and females.
        Otherwise, random pairing is performed. By default, False.
    exhaustive : bool, optional
        If True, perform exhaustive enumeration of potential mates. If False, perform random sampling. By default, True.
    """
    def __init__(self,
                 offspring_per_pair: xft.utils.VariableCount = xft.utils.ConstantCount(1),
                 mates_per_female: xft.utils.VariableCount =  xft.utils.ConstantCount(2),
                 female_offspring_per_pair: xft.mate.Union[str, xft.utils.VariableCount] = 'balanced', ## doesn't make total sense
                 sex_aware: bool = False,
                 exhaustive: bool = True,
                 ):
        super().__init__(self,
                         offspring_per_pair=offspring_per_pair,
                         mates_per_female=mates_per_female,
                         female_offspring_per_pair=female_offspring_per_pair,
                         sex_aware=sex_aware,
                         exhaustive=exhaustive,
                         component_index = None,
                         haplotypes = False,
                         )


    def mate(self,
             haplotypes: xr.DataArray = None,
             phenotypes: xr.DataArray = None,
             control: dict = None,
             ):
        """
        Mate individuals randomly with balanced numbers of males and females.

        Parameters
        ----------
        haplotypes : xr.DataArray, optional
            Array containing haplotypes, by default None
        phenotypes : xr.DataArray, optional
            Array containing phenotypes, by default None
        control : dict, optional
            Control dictionary, by default None

        Returns
        -------
        MateAssignment
            An object containing the maternal and paternal sample indices, the number of offspring per pair,
            and the number of female offspring per pair.
        """
        self._sample_indexer = haplotypes.xft.get_sample_indexer()

        if self.sex_aware:
            female_indices = haplotypes.sample[haplotypes.sex == 0]
            male_indices = haplotypes.sample[haplotypes.sex == 1]
            if len(female_indices) != len(male_indices):
                print("Warning: unequal sex ratio.")
        else:
            female_indices, male_indices = self.get_potential_mates(haplotypes, phenotypes)


        female_indices = np.random.permutation(female_indices)
        male_indices = np.random.permutation(male_indices)
        if len(female_indices)!= len(male_indices):
            # resample the smallest to get equal numbers
            if len(female_indices) < len(male_indices):
                # resample
                female_indices = np.random.choice(female_indices,len(male_indices),replace=True)
            else:
                male_indices = np.random.choice(male_indices, len(female_indices), replace=True)
            print(f"Size of m/f mating pool: {len(male_indices)} {len(female_indices)}.")


        return self.enumerate_assignment(female_indices=female_indices,
                                         male_indices=male_indices,
                                         haplotypes=haplotypes,
                                         phenotypes=phenotypes,
                                         )

def compute_pheno_covariances(
    simulation,
) -> pd.DataFrame:
    """
    Build a phenotype-only covariance table for the current generation.

    Returns a DataFrame with upper-triangle entries for:
    - pedigree_id1, pedigree_id2
    - Y: phenotype covariance
    """
    # Identify individuals in the current generation.
    ids = list(simulation.pedigree.generation(simulation.generation).nodes)

    # Extract and standardize phenotypes.
    pd_pheno = simulation.phenotypes.xft.as_pd()
    pd_pheno = pd_pheno.sort_index(axis=1)
    # normalize phenotype values
    # use .loc to avoid chained indexing; chained indexing will break in pandas 3.0.
    # iterate over all phenotype components through hierarchical indexing
    for pheno in pd_pheno.columns.get_level_values('phenotype_name').unique():
        for comp in pd_pheno.columns.get_level_values('component_name').unique():
            mask = (pd_pheno.columns.get_level_values('phenotype_name') == pheno) & \
                   (pd_pheno.columns.get_level_values('component_name') == comp)
            if mask.any():
                pd_pheno.loc[:, mask] = (pd_pheno.loc[:, mask] - pd_pheno.loc[:, mask].mean()) / pd_pheno.loc[:, mask].std()
    # Phenotype-related covariance matrices.
    y_vals = pd_pheno.loc[:, ('Y', 'phenotype')].values
    pcov = y_vals @ y_vals.T

    # Upper-triangle entries (including diagonal).
    num_individuals = len(ids)
    tri_idx = np.triu_indices(num_individuals, k=0)
    individual_pairs = [
        (ids[i], ids[j])
        for i in range(num_individuals)
        for j in range(i, num_individuals)
    ]

    data = {
        "xftsim_id1": [pair[0] for pair in individual_pairs],
        "xftsim_id2": [pair[1] for pair in individual_pairs],
        "phenotype_covariance": pcov[tri_idx],
    }

    covariances = pd.DataFrame(data)
    return covariances

if h2 + b2 > 1.0:
    raise ValueError("Require h2 + b2 <= 1.0.")

arch, _ = vct_architecture(
    founders,
    h2=h2,
    b2=b2,
    parental_coeff=parental_coef,
    n=n_indivs,
    n_causal=n_causal,
    w=0,
)

mating_regime = my_RandomMatingRegime(
    mates_per_female=1,
    offspring_per_pair=xft.utils.PoissonCount(2),
    sex_aware=True,
    female_offspring_per_pair = 'balanced'
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

# Run generations with timings
num_generations = int(snakemake.wildcards['generations'])
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Running ' + str(num_generations) + ' generations...')
sim.run(num_generations)
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Complete.')
    
# Write genotypes to plink format for GRM computation with timings
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Writing to plink format...')
arch = snakemake.wildcards.get('arch', arch)
gen = snakemake.wildcards.get('generations', num_generations)
xft.io.write_to_plink1(sim.haplotypes, f"data/{arch}/plink/h2_{h2}_b2_{b2}_pc_{parental_coef}_gen_{gen}_rep{rep}_plink")
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Complete.')

# Compute phenotypic covariances
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Computing phenotype covariances...')
pcov_df = compute_pheno_covariances(sim)
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Complete.')

# Write phenotypic covariances
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Writing phenotype covariances...')
pcov_df.to_parquet(snakemake.output['phenotype_covariances'], index=False)
now = datetime.datetime.now()
print('[' + str(now) + ']' + ' Complete.')