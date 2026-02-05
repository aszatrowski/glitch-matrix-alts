# h2_VALUES = [0.001, 0.5, 1]
# b2_VALUES = [0.0, 0.5, 1]
h2_VALUES = [0.001]
b2_VALUES = [0.5]
N_REPLICATES = 2

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
            outputs.append(f"data/vct/covariances_h2_{h2}_b2_{b2}_merged.csv")
            # outputs.append(f"figures/vct/covariances_binned_h2_{h2}_b2_{b2}_rep{rep}.png")
    return outputs

rule all:
    input: 
        get_all_outputs_vct()

rule sim_vct:
    output: 
        covariances_csv = "data/vct/covariances_h2_{h2}_b2_{b2}_rep{rep}.csv"
    params:
        n_indivs = 1000,
        m_variants = 500,
        n_causal = 100
    log:
        "logs/vct/covariances_h2_{h2}_b2_{b2}_rep{rep}.log"
    conda: "envs/shared-e-env.yaml"
    script:
        "scripts/sim_vct.py" 

rule merge_replicates:
    input:
        replicates = expand("data/vct/covariances_{{params}}_rep{rep}.csv", rep=range(N_REPLICATES))
    output:
        merged_replicates = "data/vct/covariances_{params}_merged.csv"
    conda: "envs/r-tools.yaml"
    script: "scripts/merge_replicates.R"

# rule plot_pheno_covariance_binned:
#     input: 
#         covariances_csv = "data/{arch}/covariances_h2_{h2}_b2_{b2}.csv"
#     output: 
#         covariance_plot = "figures/{arch}/covariances_binned_h2_{h2}_b2_{b2}.png",
#     conda: "envs/r-tools.yaml"
#     script: "scripts/plot_covariance_binned.R"

# rule plot_pheno_covariance_all:
#     input: 
#         covariances_csv = "data/{arch}/covariances_h2_{h2}_b2_{b2}.csv"
#     output: 
#         covariance_plot = "figures/{arch}/covariances_h2_{h2}_b2_{b2}.png",
#     conda: "envs/r-tools.yaml"
#     script: "scripts/plot_covariance_all.R"