rule all:
    input: 
        expand(
            "figures/covariance_{arch}_{h2}_all.png",
            arch = "GCTA",
            h2 = [0.5, 1],
        ),
        expand(
            "figures/covariance_{arch}_{h2}_binned.png",
            arch = "GCTA",
            h2 = [0.5, 1]
        )

rule sim:
    output: 
        covariances_csv = "data/covariances_{arch}_{h2}.csv"
    params:
        n_indivs = 1000,
        m_variants = 1000
    conda: "envs/shared-e-env.yaml"
    script: "scripts/sim.py"

rule plot_pheno_covariance_binned:
    input: 
        covariances_csv = "data/covariances_{arch}_{h2}.csv"
    output: 
        covariance_plot = "figures/covariance_{arch}_{h2}_binned.png"
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_all:
    input: 
        covariances_csv = "data/covariances_{arch}_{h2}.csv"
    output: 
        covariance_plot = "figures/covariance_{arch}_{h2}_all.png",
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_all.R"