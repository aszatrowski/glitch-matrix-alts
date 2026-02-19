rule sim_vct:
    output: 
        covariances_plink = temp(multiext(
            "data/{arch}/plink/h2_{h2}_b2_{b2}_pc_{parental_coef}_rep{rep}_plink",
            ".bed", ".bim", ".fam"
        )),
        phenotype_covariances = temp("data/{arch}/pcov/h2_{h2}_b2_{b2}_pc_{parental_coef}_rep{rep}.parquet")
    params:
        n_indivs = 1500,
        m_variants = 5e4,
        n_causal = 5e4 - 1,
        chrom_count = 2,
        generations = 10
    log:
        "logs/{arch}/h2_{h2}_b2_{b2}_pc_{parental_coef}_rep{rep}.log"
    resources:
        mem = "48G",
        runtime = 15
    conda: "envs/xftsim.yaml"
    script: "scripts/sim_vct.py" 

rule plink_compute_grm:
    input: 
        covariances_plink = multiext(
            "data/{arch}/plink/h2_{h2}_b2_{b2}_pc_{parental_coef}_rep{rep}_plink",
            ".bed", ".bim", ".fam"
        )
    output: 
        temp(multiext(
            "data/{arch}/grm/h2_{h2}_b2_{b2}_pc_{parental_coef}_rep{rep}",
            ".rel", ".rel.id", ".log"
        ))
    resources:
        mem = "8G",
        runtime = 5 
    threads: 2
    params:
        min_af = 0.01
    conda: "envs/r-plink.yaml"
    shell: 
        """
        plink2 \
        --bfile data/{wildcards.arch}/plink/h2_{wildcards.h2}_b2_{wildcards.b2}_pc_{wildcards.parental_coef}_rep{wildcards.rep}_plink \
        --min-af {params.min_af} \
        --make-rel square0 \
        --threads {threads} \
        --nonfounders \
        --out data/{wildcards.arch}/grm/h2_{wildcards.h2}_b2_{wildcards.b2}_pc_{wildcards.parental_coef}_rep{wildcards.rep}
        """

rule merge_replicates:
    input:
        grm_replicates = expand("data/{{arch}}/grm/h2_{{h2}}_b2_{{b2}}_pc_{{parental_coef}}_rep{rep}.rel", rep=range(N_REPLICATES)),
        grm_replicate_iids = expand("data/{{arch}}/grm/h2_{{h2}}_b2_{{b2}}_pc_{{parental_coef}}_rep{rep}.rel.id", rep=range(N_REPLICATES)),
        phenotype_covariance_replicates = expand("data/{{arch}}/pcov/h2_{{h2}}_b2_{{b2}}_pc_{{parental_coef}}_rep{rep}.parquet", rep=range(N_REPLICATES))
    output:
        merged_replicates = "data/{arch}/h2_{h2}_b2_{b2}_pc_{parental_coef}_covmatrix_merged.parquet"
    params:
        generations = 10
    resources:
        mem = "32G",
        runtime = 5 
    threads: 2
    conda: "envs/r-plink.yaml"
    script: "scripts/merge_replicates.R"

rule plot_pheno_covariance_binned:
    input: 
        covariances_csv = "data/{arch}/h2_{h2}_b2_{b2}_pc_{parental_coef}_covmatrix_merged.parquet"
    output: 
        covariance_plot = "figures/{arch}/h2_{h2}_b2_{b2}_pc_{parental_coef}_binned.png",
    params:
        binwidth = 0.01,
        min_obs_in_bin = 20
    resources:
        mem = "16G",
        runtime = 5 
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_all:
    input: 
        covariances_csv = "data/{arch}/h2_{h2}_b2_{b2}_pc_{parental_coef}_covmatrix_merged.parquet"
    output: 
        covariance_plot = "figures/{arch}/h2_{h2}_b2_{b2}_pc_{parental_coef}_all.png",
    resources:
        mem = "16G",
        runtime = 10 
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_all.R"

rule plot_pheno_covariance_binned_overlay:
    input: 
        covariances_csv_list = lambda wildcards: [
            f"data/{wildcards.arch}/h2_{wildcards.h2}_b2_{b2}_pc_{wildcards.parental_coef}_covmatrix_merged.parquet"
            for b2 in b2_VALUES
            if float(wildcards.h2) + float(b2) <= 1.0
        ]
    output: 
        covariance_plot = "figures/{arch}/overlay_h2_{h2}_pc_{parental_coef}_binned.png",
    params:
        binwidth = 0.01,
        min_obs_in_bin = 5
    resources:
        mem = "64G",
        runtime = 5 
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_overlay_binned.R"