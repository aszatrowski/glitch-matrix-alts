# Define your grid
h2_VALUES = [0.001, 0.25, 0.5, 1]
b2_VALUES = [0.0, 0.25, 0.5, 1]

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
        outputs.append(f"figures/vct/covariances_h2_{h2}_b2_{b2}.png")
    return outputs

rule all:
    input: 
        expand(
            "data/gcta/covariances_{h2}.csv",
            h2 = [0.5, 1],
        ),
        get_all_outputs_vct()

rule sim_gcta:
    output: 
        covariances_csv = "data/gcta/covariances_{h2}.csv"
    params:
        n_indivs = 1000,
        m_variants = 1000
    conda: "envs/shared-e-env.yaml"
    script: "scripts/sim.py"

rule sim_vct:
    output: 
        covariances_csv = "data/vct/covariances_h2_{h2}_b2_{b2}.csv"
    params:
        n_indivs = 1000,
        m_variants = 500,
        n_causal = 100
    log:
        "logs/vct/covariances_h2_{h2}_b2_{b2}.log"
    conda: "envs/shared-e-env.yaml"
    script:
        "scripts/sim_vct.py" 

# rule plot_pheno_covariance_gcta_binned:
#     input: 
#         covariances_csv = "data/covariances_gcta_{h2}.csv"
#     output: 
#         covariance_plot = "figures/covariance_gcta_{h2}_binned.png"
#     conda: "envs/r-tools.yaml"
#     script: "scripts/plot_covariance_binned.R"

# rule plot_pheno_covariance_gcta_all:
#     input: 
#         covariances_csv = "data/covariances_gcta_{h2}.csv"
#     output: 
#         covariance_plot = "figures/covariance_gcta_{h2}_all.png",
#     conda: "envs/r-tools.yaml"
#     script: "scripts/plot_covariance_all.R"

rule plot_pheno_covariance_vct_binned:
    input: 
        covariances_csv = "data/vct/covariances_h2_{h2}_b2_{b2}.csv"
    output: 
        covariance_plot = "figures/vct/covariances_binned_h2_{h2}_b2_{b2}.png",
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_vct_all:
    input: 
        covariances_csv = "data/vct/covariances_h2_{h2}_b2_{b2}.csv"
    output: 
        covariance_plot = "figures/vct/covariances_h2_{h2}_b2_{b2}.png",
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_all.R"