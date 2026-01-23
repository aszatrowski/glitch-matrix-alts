rule all:
    input: 
        # expand(
        #     "results/covariances_{arch}.png",
        #     arch = "GCTA"
        # )
        expand(
            "data/covariances_{arch}.csv",
            arch = "GCTA"
        )

rule sim:
    output: 
        "data/covariances_{arch}.csv"
    conda: "envs/shared-e-env.yaml"
    script: "scripts/sim.py"