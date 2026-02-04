# Data assembly and statistical analysis for HTT project

Having extracted a list of clinical terms from the UKBB (available in: `dating_and_threshold_ukbb_columns.tsv`), we initially assembled the genotype data and clinical data
and identified ICD10 codes that were not formally dated and were present >= 2 times in people with >=36 CAG repeats at the HTT locus. This was done using:

```
Rscript extracting_icd10_codes_for_dating.R #outputs icds_gt36.tsv and ukbb_informal_dated_main_icd10_codes.tsv
```

The resulting files were converted into events using python:
```
python parse_icd10_dates.py 
```
When informally dated ICD10 codes had been assembled into a file- the assembly and dating of all genotype data, formal neurological diagnoses, relvant cognitive features
and relevant self-reported data were assembled using the following scripts:

```
Rscript bind_primary_clinical_terms_to_gt_data.R
Rscript HTT_study_cog_scores.R
Rscript HTT_date_self_reported.R
```
Another set of fields was requested at this point in the project (listed in `ukbb_encodings.txt`)- these were assembled using:
```
Rscript assemble_new_columns.R
```

Next a first-pass statistical analysis of events enriched in different cohorts was done using:
```
Rscript HTT_assessing_rates_across_stages.R
```
