h2_VALUES = ["0.0001", "1.0"]
b2_VALUES = ["0.0", "0.25", "0.5", "0.75", "0.9999"]
N_REPLICATES = 50

def get_valid_combinations():
    valid = []
    for h2 in h2_VALUES:
        for b2 in b2_VALUES:
            if float(h2) + float(b2) <= 1.0:  # Coerce for comparison
                valid.append((h2, b2))  # But keep as strings
    return valid

def get_all_outputs_vct():
    outputs = []
    for h2, b2 in get_valid_combinations():
        for rep in range(N_REPLICATES):
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_all.png")
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_binned.png")
    return outputs

def get_overlay_outputs():
    outputs = []
    for h2, b2 in get_valid_combinations():
        outputs.append(f"figures/vct/overlay_h2_{h2}_b2_binned.png")
    return outputs

rule all:
    input: 
        get_all_outputs_vct(),
        get_overlay_outputs()
    
rule sim_vct:
    output: 
        covariances_plink = temp(multiext(
            "data/vct/plink/h2_{h2}_b2_{b2}_rep{rep}_plink",
            ".bed", ".bim", ".fam"
        )),
        phenotype_covariances = temp("data/vct/pcov/h2_{h2}_b2_{b2}_rep{rep}.parquet")
    params:
        n_indivs = 1500,
        m_variants = 5e5,
        n_causal = 5e5 - 1,
        chrom_count = 2,
        generations = 10
    log:
        "logs/vct/h2_{h2}_b2_{b2}_rep{rep}.log"
    resources:
        mem = "48G",
        runtime = 15
    conda: "envs/xftsim.yaml"
    script: "scripts/sim_vct.py" 

rule plink_compute_grm:
    input: 
        covariances_plink = multiext(
            "data/vct/plink/h2_{h2}_b2_{b2}_rep{rep}_plink",
            ".bed", ".bim", ".fam"
        )
    output: 
        temp(multiext(
            "data/vct/grm/h2_{h2}_b2_{b2}_rep{rep}",
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
        --bfile data/vct/plink/h2_{wildcards.h2}_b2_{wildcards.b2}_rep{wildcards.rep}_plink \
        --min-af {params.min_af} \
        --make-rel square0 \
        --threads {threads} \
        --nonfounders \
        --out data/vct/grm/h2_{wildcards.h2}_b2_{wildcards.b2}_rep{wildcards.rep}
        """

rule merge_replicates:
    input:
        grm_replicates = expand("data/vct/grm/h2_{{h2}}_b2_{{b2}}_rep{rep}.rel", rep=range(N_REPLICATES)),
        grm_replicate_iids = expand("data/vct/grm/h2_{{h2}}_b2_{{b2}}_rep{rep}.rel.id", rep=range(N_REPLICATES)),
        phenotype_covariance_replicates = expand("data/vct/pcov/h2_{{h2}}_b2_{{b2}}_rep{rep}.parquet", rep=range(N_REPLICATES))
    output:
        merged_replicates = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.parquet"
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
        covariances_csv = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.parquet"
    output: 
        covariance_plot = "figures/vct/h2_{h2}_b2_{b2}_binned.png",
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
        covariances_csv = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.parquet"
    output: 
        covariance_plot = "figures/vct/h2_{h2}_b2_{b2}_all.png",
    resources:
        mem = "16G",
        runtime = 10 
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_all.R"

rule plot_pheno_covariance_binned_overlay:
    input: 
        covariances_csv_list = lambda wildcards: [
            f"data/vct/h2_{wildcards.h2}_b2_{b2}_covmatrix_merged.parquet"
            for b2 in b2_VALUES
            if float(wildcards.h2) + float(b2) <= 1.0
        ]
    output: 
        covariance_plot = "figures/vct/overlay_h2_{h2}_b2_binned.png",
    params:
        binwidth = 0.01,
        min_obs_in_bin = 5
    resources:
        mem = "64G",
        runtime = 5 
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_overlay_binned.R"