def get_cellranger(wildcards):
    return os.path.join(sampledir, samples.loc[wildcards.sample, "cellranger"])

def get_rawh5(wildcards):
    x = samples.loc[wildcards.sample, "cellranger"]
    return os.path.join(sampledir, f"{x}/raw_feature_bc_matrix.h5")

def get_cellrangersummary(wildcards):
    x = samples.loc[wildcards.sample, "cellranger"]
    return (wildcards.sample, os.path.join(sampledir, f"{x}/metrics_summary.csv"))
