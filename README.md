# PPARGC1B kidney state analysis

Reproducibility materials for **PPARGC1B reduction and heterogeneous metabolic programs across human proximal tubule states**.

## Scope
Three related human state contrasts are analysed separately; no combined effect or treatment/progression inference is made. GEO GSE254185, GSE195460 and the deposited GSE183277 snCv3 matrix are the human sources. GSE139107 supplies descriptive mouse support. Shared donor identity across studies cannot be verified.

## Run the core analysis
Use Python with numpy, pandas and scipy, and R with edgeR, limma and statmod. From this repository root:

```sh
python core/analyze_programs.py
python core/analyze_external.py
python core/analyze_segments.py
Rscript core/fit_edger.R analysis_inputs
```

External edgeR uses its script's working directory. Run `Rscript ../core/fit_external.R` from `analysis_inputs`. All-gene PT donor pseudobulks and count-model inputs are supplied. These are derived from public matrices; nuclei are not biological replicates. The validation compressed table includes only PT and PTVCAM1 states required by the runnable core; all-feature UMI denominators are preserved.

## Provenance and limits
`analysis_inputs/scripts` preserves historical scripts, including plotting and upstream downloads. These scripts retain their original project-relative layout and may require path adaptation and downloading the public source matrices. They are not all invoked by the core commands. The compact archive supplies exact gene lists, source URLs, count matrices, numerical outputs, eligibility records and R session logs. It does not redistribute the large original sequencing matrices. Source publications and third-party data retain their own terms; this bundle does not relicense their data. AI assistance covered analysis-script development/execution and manuscript preparation, and is disclosed in the manuscript.

No new human or animal experiments were performed. Subgroup and segment sensitivities are exploratory; CKD PT-S2 is nonsignificant.

## GitHub distribution

Download GitHub_code_and_derived_inputs.zip from this repository and extract it into a local working directory. It contains the complete core/ and analysis_inputs/ directories, including derived public-source inputs. The core scripts are also available as browsable files in core/; the compressed bundle supplies their required matrices. Run the commands above from the extracted directory.

The final submission figure numbering and exact source tables are indexed in figure_source_index.csv. Main Figure 4 is external validation, main Figure 5 is exploratory context, and ATAC intervals are Supplementary Figure S3. No original sequencing matrices or participant-identifying records are uploaded.
