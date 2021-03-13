def get_cellranger(wildcards):
    return samples.loc[wildcards.sample, "cellranger"]


def get_rawh5(wildcards):
    x = samples.loc[wildcards.sample, "cellranger"]
    return f"{x}/raw_feature_bc_matrix.h5"
