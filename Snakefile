from src.distribution_sep import AssociationGroupAdjustment
import os
# CONFIG
test = ["as2.assoc", ]


# FUNCTIONS
def input_plots(wildcards):
    checkpoint_output = checkpoints.generate_plots.get(**wildcards).output[0]
    print(checkpoint_output)
    plots = os.listdir(checkpoint_output)
    plots = [p for p in plots if p[-4:-1] == ".pn"]
    return expand(checkpoint_output + "/{image}",
            image = plots)


# RULE
rule all:
    input:
        expand("results/{assoc}.html", assoc=test)


rule generate_infographics:
    input:
        input_plots
    output:
        html_report = "results/{assoc}.html"
    shell:
        """
        mkdir -p results
        Rscript -e "library(rmdformats); rmarkdown::render('src/generate_sample_report.Rmd', output_file='../{output.html_report}', output_format=c('readthedown'))" \
            --args {wildcards.assoc} {input}
        """


checkpoint generate_plots:
    input:
        feature_matrix = "run_folder/{assoc}.mat",
        annotated      = "run_folder/{assoc}.ann"
    output:
        plot_dir = directory("run_folder/plots/{assoc}")
    run:
        shell(f"mkdir -p run_folder/plots/{wildcards.assoc}")
        as_obj = AssociationGroupAdjustment(
            input.annotated,
            input.feature_matrix,
            f"run_folder/plots/{wildcards.assoc}/")
        res = as_obj.get_fdr_per_group(True)


rule get_feature_matrix:
    input:
        "input/{assoc}.tsv"
    output:
        feature_matrix = "run_folder/{assoc}.mat",
        annotated = "run_folder/{assoc}.annotated.tsv"
    shell:
        """
        Rscript src/get_gene_from_rsid.R \
            --args
            --association_file={input} \
            --annotated_output={output.annotated} \
            --feature_matrix_output={output.feature_matrix}
        """
