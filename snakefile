h2_VALUES = ["0.0001", "1.0"]
b2_VALUES = ["0.0", "0.25", "0.5", "0.75", "0.9999"]
parental_coef_VALUES = ["0.5"]

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
        for parental_coef in parental_coef_VALUES:
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_pc_{parental_coef}_all.png")
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_pc_{parental_coef}_binned.png")
    return outputs

def get_overlay_outputs():
    outputs = []
    for h2, b2 in get_valid_combinations():
        for parental_coef in parental_coef_VALUES:
            outputs.append(f"figures/vct/overlay_h2_{h2}_pc_{parental_coef}_binned.png")
    return outputs

include: "workflow/vct.smk"

rule all:
    input: 
        get_all_outputs_vct(),
        get_overlay_outputs()
    