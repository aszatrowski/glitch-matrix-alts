rule all:
    input: 
        expand(
            "figures/covariance_{arch}.png",
            arch = "GCTA"
        ),
        expand(
            "data/covariances_{arch}.csv",
            arch = "GCTA"
        )

rule sim:
    output: 
        covariances_csv = "data/covariances_{arch}.csv"
    conda: "envs/shared-e-env.yaml"
    script: "scripts/sim.py"

rule plot_pheno_covariance:
    input: 
        covariances_csv = "data/covariances_{arch}.csv"
    output: 
        covariance_plot = "figures/covariance_{arch}.png",
    conda: "envs/shared-e-env.yaml"
    script: "scripts/plot_covariance.R"