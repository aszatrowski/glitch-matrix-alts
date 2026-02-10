h2_VALUES = ["0.001", "0.5", "1.0"]
b2_VALUES = ["0.0", "0.25", "0.5", "0.75", "1.0"]
N_REPLICATES = 2

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

rule all:
    input: 
        get_all_outputs_vct()

rule sim_vct:
    output: 
        covariances_plink = temp(multiext(
            "data/vct/plink/h2_{h2}_b2_{b2}_rep{rep}_plink",
            ".bed", ".bim", ".fam"
        )),
        phenotype_covariances = temp("data/vct/pcov/h2_{h2}_b2_{b2}_rep{rep}.csv")
    params:
        n_indivs = 1200,
        m_variants = 300,
        n_causal = 200,
        generations = 10
    log:
        "logs/vct/h2_{h2}_b2_{b2}_rep{rep}.log"
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
            ".rel", ".rel.id"
        ))
    log:
        "logs/vct/plink/h2_{h2}_b2_{b2}_rep{rep}.log"
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
        --out data/vct/grm/h2_{wildcards.h2}_b2_{wildcards.b2}_rep{wildcards.rep} \
        &> {log}
        """

rule merge_replicates:
    input:
        grm_replicates = expand("data/vct/grm/h2_{{h2}}_b2_{{b2}}_rep{rep}.rel", rep=range(N_REPLICATES)),
        grm_replicate_iids = expand("data/vct/grm/h2_{{h2}}_b2_{{b2}}_rep{rep}.rel.id", rep=range(N_REPLICATES)),
        phenotype_covariance_replicates = expand("data/vct/pcov/h2_{{h2}}_b2_{{b2}}_rep{rep}.csv", rep=range(N_REPLICATES))
    output:
        merged_replicates = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.csv"
    params:
        generations = 10
    threads: 2
    conda: "envs/r-plink.yaml"
    script: "scripts/merge_replicates.R"

rule plot_pheno_covariance_binned:
    input: 
        covariances_csv = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.csv"
    output: 
        covariance_plot = "figures/vct/h2_{h2}_b2_{b2}_binned.png",
    params:
        binwidth = 0.01,
        min_obs_in_bin = 5
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_binned.R"

rule plot_pheno_covariance_all:
    input: 
        covariances_csv = "data/vct/h2_{h2}_b2_{b2}_covmatrix_merged.csv"
    output: 
        covariance_plot = "figures/vct/h2_{h2}_b2_{b2}_all.png",
    conda: "envs/r-plink.yaml"
    script: "scripts/plot_covariance_all.R"
