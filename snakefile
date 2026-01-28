rule all:
    input: 
        expand(
            "figures/covariance_{arch}_{h2}_{binning}.png",
            arch = "GCTA",
            h2 = [0.5, 1],
            binning = ['all', 'binned']
        )
        # ),
        # expand(
        #     "data/covariances_{arch}_{h2}.csv",
        #     arch = "GCTA"
        # )

rule sim:
    output: 
        covariances_csv = "data/covariances_{arch}_{h2}.csv"
    conda: "envs/shared-e-env.yaml"
    script: "scripts/sim.py"

rule plot_pheno_covariance_binned:
    input: 
        covariances_csv = "data/covariances_{arch}_{h2}.csv"
    output: 
        covariance_plot = "figures/covariance_{arch}_{0.5}_binned.png",
    conda: "envs/shared-e-env.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_all:
    input: 
        covariances_csv = "data/covariances_{arch}_{h2}.csv"
    output: 
        covariance_plot = "figures/covariance_{arch}_{h2}_{binning}.png",
    conda: "envs/shared-e-env.yaml"
    script: "scripts/plot_covariance_all.R"