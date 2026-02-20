h2_VALUES = ["0.0001", "0.5", "1.0"]
b2_VALUES = ["0.0", "0.25", "0.5", "0.75", "0.9999"]
parental_coef_VALUES = ["0.0", "0.25", "0.5", "0.75", "1.0"]

N_REPLICATES = 10
N_VARIANTS = 1e5
N_CAUSAL_VARIANTS = N_VARIANTS - 1

def get_valid_combinations():
    valid = []
    for h2 in h2_VALUES:
        for b2 in b2_VALUES:
            if float(h2) + float(b2) <= 1.0:
                valid.append((h2, b2))
    return valid

def get_valid_parental_coefs(b2):
    """
    If b2 is 0, changing parental_coef will make no difference, so only run once with parental_coef = 0.5.
    """
    if float(b2) == 0.0:
        return ["0.5"]
    return parental_coef_VALUES

def get_all_outputs_vct():
    outputs = []
    for h2, b2 in get_valid_combinations():
        for parental_coef in get_valid_parental_coefs(b2):
            outputs.append(f"figures/vct/h2_{h2}_b2_{b2}_pc_{parental_coef}_all.png")
    return outputs

def get_overlay_outputs():
    outputs = []
    for h2, b2 in get_valid_combinations():
        for parental_coef in get_valid_parental_coefs(b2):
            outputs.append(f"figures/vct/b2_overlay/overlay_h2_{h2}_pc_{parental_coef}_binned.png")
    return outputs

include: "workflow/vct.smk"

rule all:
    input: 
        get_overlay_outputs(),
        get_all_outputs_vct()