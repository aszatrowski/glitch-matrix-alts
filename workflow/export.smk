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
        "figures/vct/b2_overlay/overlay_h2_0.0001_pc_1.0_gen_5_binned.png",
        "figures/vct/b2_overlay/overlay_h2_0.0001_pc_1.0_gen_10_binned.png",
        "figures/vct/b2_overlay/overlay_h2_0.0001_pc_1.0_gen_20_binned.png",
    output: 
        pptx = "export/results.pptx"
    conda:
        "envs/render.yaml"
    script: 
        "scripts/make_presentation.py"