rule parental_imbalance_theory:
    input:
        rmd = "export/theory.Rmd"
    output: 
        html = "export/theory.pdf"
    conda:
        "envs/render.yaml"
    shell:
        """
        Rscript -e \"if (!tinytex::is_tinytex()) tinytex::install_tinytex()\"
        Rscript -e \"rmarkdown::render('{input}', output_file='$(basename {output})', knit_root_dir=getwd(), clean=TRUE)\"
        """

rule pptx:
    input: 
        "figures/vct/h2_0.0001_b2_0.0_pc_0.5_gen_2_all.png",
        "figures/vct/h2_0.0001_b2_0.0_pc_0.5_gen_10_all.png",
    output: 
        "export/results.pptx"
    conda:
        "envs/render.yaml"
    script: 
        "scripts/make_presentation.py"