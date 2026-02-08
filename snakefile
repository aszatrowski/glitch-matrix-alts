h2_VALUES = [0.001, 0.5, 1]
b2_VALUES = [0.0, 0.25, 0.5, 0.75, 1]
N_REPLICATES = 50

def get_valid_combinations():
    valid = []
    for h2 in h2_VALUES:
        for b2 in b2_VALUES:
            if h2 + b2 <= 1.0:  # Ensure environmental component is non-negative
                valid.append((h2, b2))
    return valid

def get_all_outputs_vct():
    outputs = []
    for h2, b2 in get_valid_combinations():
        for rep in range(N_REPLICATES):
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_all.png")
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_binned.png")
    return outputs

rule all:
    input: 
        get_all_outputs_vct()

rule sim_vct:
    output: 
        covariances_csv = temp("data/vct/h2_{h2}_b2_{b2}_rep{rep}_covmatrix.csv")
    params:
        n_indivs = 1200,
        m_variants = 300,
        n_causal = 200,
        generations = 10
    log:
        "logs/vct/h2_{h2}_b2_{b2}_rep{rep}.log"
    conda: "envs/xftsim.yaml"
    script: "scripts/sim_vct.py" 

rule merge_replicates:
    input:
        replicates = expand("data/vct/h2_{{h2}}_b2_{{b2}}_rep{rep}_covmatrix.csv", rep=range(N_REPLICATES))
    output:
        merged_replicates = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.csv"
    conda: "envs/r-tools.yaml"
    script: "scripts/merge_replicates.R"

rule plot_pheno_covariance_binned:
    input: 
        covariances_csv = "data/{arch}/h2_{h2}_b2_{b2}_covmatrix_merged.csv"
    output: 
        covariance_plot = "figures/{arch}/h2_{h2}_b2_{b2}_binned.png",
    params:
        binwidth = 0.01,
        min_obs_in_bin = 5
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_all:
    input: 
        covariances_csv = "data/{arch}/h2_{h2}_b2_{b2}_covmatrix_merged.csv"
    output: 
        covariance_plot = "figures/{arch}/h2_{h2}_b2_{b2}_all.png",
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_all.R"
