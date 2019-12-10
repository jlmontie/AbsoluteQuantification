# Synergy Validation
This directory contains work for the Synergy validation. Quantification
was performed on simulated reads from 16S/18S of these organisms:
- *Escherichia coli*
- *Klebsiella pneumoniae*
- *Candida albicans*
- *Candida tropicalis*

Validation was performed with the following steps:
1. Gather reference sequences from our internal database. See `ref_seqs`
directory. The references were gathered with
`get_selected_target_ref_seqs.py` on the CHPC rubi server. Absolute paths
are provided in the script.

2. Reads were simulated from the reference sequences. See `simulated_reads`
directory. Reads were simulated with the `simulate_reads.py` script which
calls simulation scripts in the `taxonomer` repo. Reads were simulated at
depths indicated in the validation plan.

3. The simulated reads were combined into a single fasta in the proportions
designated in the validation plan. See `simulated_reads` directory,
`combine_reads.py`.

4. Simulated reads were processed with the `synergy_classification` pipeline.
See `pipeline` directory. The input file was prepared with
`write_input_yaml.py`. The pipeline was called with
`run_simulated_reads_classification_pipeline.sh`.

5. The results were analyzed. See `report` directory. The pipeline output was
analyzed with `analyze_in_silico_validation_run.py`. This produced
`in_silico_quantification_results.csv` which was manually adjusted for
compatibility with Latex for the final report
(`in_silico_quantification_final_report_superscript.csv`). This directory also
containse a script for generating an explanatory figure for the validation
report describing the quantification algorithm.

The directories `ge_cfu_correlation` and `plan` include exploratory analysis
for the validation plan. `Pipfile` in the parent `AbsoluteQuantification`
directory is a requirements file produced by `pipenv` and provides the
packages required to run these scripts. Run `pipenv shell` to activate the
pipenv environment.