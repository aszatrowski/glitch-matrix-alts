rule all:
    input: 
        expand(
            "data/gcta/covariances_{h2}.csv",
            h2 = [0.5, 1],
        ),
        expand(
            "figures/covariance_vct_{h2}_{b2}.png",
            h2 = [0.001, 0.5, 1],
            b2 = [0.1, 0.25]
        ),
        # expand(
        #     "figures/covariance_vct_{h2}_{b2}_binned.png",
        #     h2 = [0.5, 1],
        #     b2 = [0, 0.25]
        # )

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
        covariances_csv = "data/vct/covariances_{h2}_{b2}.csv"
    params:
        n_indivs = 1000,
        m_variants = 1000,
        n_causal = 50
    conda: "envs/shared-e-env.yaml"
    script: "scripts/sim.py"

rule plot_pheno_covariance_gcta_binned:
    input: 
        covariances_csv = "data/covariances_gcta_{h2}.csv"
    output: 
        covariance_plot = "figures/covariance_gcta_{h2}_binned.png"
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_gcta_all:
    input: 
        covariances_csv = "data/covariances_gcta_{h2}.csv"
    output: 
        covariance_plot = "figures/covariance_gcta_{h2}_all.png",
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_all.R"

rule plot_pheno_covariance_vct_binned:
    input: 
        covariances_csv = "data/vct/covariances_{h2}_{b2}.csv"
    output: 
        covariance_plot = "figures/covariance_vct_{h2}_{b2}_binned.png"
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_vct_all:
    input: 
        covariances_csv = "data/vct/covariances_{h2}_{b2}.csv"
    output: 
        covariance_plot = "figures/covariance_vct_{h2}_{b2}.png",
    conda: "envs/r-tools.yaml"
    script: "scripts/plot_covariance_all.R"