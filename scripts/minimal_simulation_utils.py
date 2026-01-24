"""
Minimal helpers for the VCT demo.

This module intentionally contains only what the minimal simulation needs:
GRM construction and a covariance table builder for the final generation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_grm(simulation, maf_min: float = 0.01) -> np.ndarray:
    """
    Compute a genetic relationship matrix (GRM) from a simulation object.

    Assumes the simulation exposes:
    - simulation.current_afs_empirical
    - simulation.current_std_genotypes (shape: n_individuals x n_variants)
    """
    afs = simulation.current_afs_empirical
    mafs = np.minimum(afs, 1 - afs)
    maf_mask = (mafs >= maf_min) & (mafs <= 0.5)

    geno = simulation.current_std_genotypes[:, maf_mask]
    var = np.nanvar(geno, axis=0)
    var_mask = (var > 0) & ~np.isnan(var)
    geno = geno[:, var_mask]

    return np.cov(geno)


def get_pedigree_distance(sim, maximum_depth=None) -> pd.DataFrame:
    """
    Trace pedigree relationships back and return pairwise distances.

    Uses sgkit to compute additive relationship (kinship) on a pedigree
    derived from the simulation's pedigree graph.
    """
    import networkx as nx
    import xarray as xr
    import sgkit as sg

    pedigree = sim.pedigree.G

    # Convert the graph into parent-child trios (sgkit format).
    trios = list(nx.compute_v_structures(pedigree))
    trio_df = pd.DataFrame(trios, columns=["parent1", "child", "parent2"])
    top_gens = list(nx.topological_generations(pedigree))
    ng = len(top_gens)

    if maximum_depth is None or maximum_depth >= ng:
        maximum_depth = ng
        founders = top_gens[0]
        founder_df = pd.DataFrame(founders, columns=["child"])
        founder_df["parent1"] = "."
        founder_df["parent2"] = "."
        trio_df = pd.concat([trio_df, founder_df], ignore_index=True)
        print("Searching the full depth pedigree.")
    else:
        gens = top_gens[-maximum_depth:]
        trio_df = trio_df[trio_df["child"].isin([item for sub in gens for item in sub])]
        founders = top_gens[-(maximum_depth + 1)]
        founder_df = pd.DataFrame(founders, columns=["child"])
        founder_df["parent1"] = "."
        founder_df["parent2"] = "."
        trio_df = pd.concat([trio_df, founder_df], ignore_index=True)
        print("Searching partial depth pedigree.")

    pedigree_sg = trio_df.astype("U").set_index("child", drop=True)
    parent_id = xr.DataArray(pedigree_sg, dims=("samples", "parents"))
    ped = xr.Dataset()
    ped["sample_id"] = parent_id.samples
    ped["parent_id"] = parent_id
    ped = sg.pedigree_kinship(
        ped,
        allow_half_founders=True,
        return_relationship=True,
    )
    ped = ped.assign_coords(dict(samples_0=ped.samples.values, samples_1=ped.samples.values))
    ped = ped.compute()

    # Use the final generation (leaf) samples for pairwise distances.
    leaves = top_gens[-2]
    ped_leaves = ped.stat_pedigree_relationship.sel(samples_0=leaves, samples_1=leaves)
    ped_df = ped_leaves.to_dataframe().reset_index()
    ped_df.columns = ["pedigree_id1", "pedigree_id2", "Pedigree_Distance"]
    return ped_df


def make_covariance_matrix(
    simulation,
    maf: float = 0.01,
    include_pedigree: bool = True,
    pedigree_depth=None,
) -> pd.DataFrame:
    """
    Build a covariance table for the current generation.

    Returns a DataFrame with upper-triangle entries for:
    - pedigree_id1, pedigree_id2
    - Y: phenotype covariance
    - A: additive genetic covariance
    - C: heritable noise covariance
    - Genetic_Covariance: GRM entry
    - Pedigree_Distance (optional, from sgkit)
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
            col = (pheno, comp, slice(None))
            pd_pheno.loc[:, col] = (pd_pheno.loc[:, col] - pd_pheno.loc[:, col].mean()) / pd_pheno.loc[:, col].std()
    # GRM for genetic covariance.
    grm = compute_grm(simulation, maf_min=maf)

    # Phenotype-related covariance matrices.
    pcov = pd_pheno.loc[:, ('Y', 'phenotype')] @ pd_pheno.loc[:, ('Y', 'phenotype')].T
    acov = pd_pheno.loc[:, ('Y', 'additiveGenetic')] @ pd_pheno.loc[:, ('Y', 'additiveGenetic')].T
    # ccov = pd.DataFrame(pd_pheno.loc[:, ('Y', 'heritableNoise')]) @ pd.DataFrame( pd_pheno.loc[:, ('Y', 'heritableNoise')]).T

    # Upper-triangle entries (including diagonal).
    num_individuals = len(ids)
    tri_idx = np.triu_indices(num_individuals, k=0)
    individual_pairs = [
        (ids[i], ids[j])
        for i in range(num_individuals)
        for j in range(i, num_individuals)
    ]

    data = {
        "pedigree_id1": [pair[0] for pair in individual_pairs],
        "pedigree_id2": [pair[1] for pair in individual_pairs],
        "Y": pcov.values[tri_idx],
        "A": acov.values[tri_idx],
        # "C": ccov.values[tri_idx],
        "Genetic_Covariance": grm[tri_idx],
    }

    covariances = pd.DataFrame(data)

    if include_pedigree:
        ped_df = get_pedigree_distance(simulation, maximum_depth=pedigree_depth)
        covariances = covariances.merge(
            ped_df,
            on=["pedigree_id1", "pedigree_id2"],
            how="left",
        )
    return covariances