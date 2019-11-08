from src.distribution_sep import AssociationGroupAdjustment
import os
# CONFIG
test = ["as2.assoc", ]


# FUNCTIONS
def input_plots(wildcards):
    checkpoint_output = checkpoints.generate_plots.get(**wildcards).output[0]
    print(checkpoint_output)
    return expand(checkpoint_output + "/{image}",
            image = os.listdir(checkpoint_output))


# RULE
rule all:
    input:
        expand("run_folder/report/{assoc}.html", assoc=test)


rule generate_infographics:
    input:
        input_plots
    output:
        html_report = "run_folder/report/{assoc}.html"
    shell:
        """
        echo {output.html_report}
        mkdir -p run_folder/report
        Rscript -e "rmarkdown::render('bin/generate_sample_report.Rmd', output_file='../{output.html_report}')" \
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
        Rscript bin/get_gene_from_rsid.R \
            --args
            --association_file={input} \
            --annotated_output={output.annotated} \
            --feature_matrix_output={output.feature_matrix}
        """