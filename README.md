# TNBC recurrence associated ceRNA code and processed data

Companion materials for **Integrative Bioinformatics Analysis Suggests a Candidate Recurrence-Associated ceRNA Network in Triple-Negative Breast Cancer**.

This release preparation contains the supplied processed miRNA results, 59 consensus miRNA targets, final GSE43502 11-gene statistics, and final ceRNA network nodes and edges. It also preserves selected original R analysis sources. The supplied final 11-gene statistics match the corresponding rows of the earlier 14-gene axis within numeric export precision.

## Scope and status

The processed final-network tables were checked successfully by the included Python validation program. Original R sources are accumulated interactive work and are not an independently tested end-to-end pipeline. R was unavailable during preparation. This package does not contain all inputs and outputs necessary to reproduce every manuscript analysis. See `documentation/RELEASE_COMPLETENESS.md` before representing this as the complete manuscript code release.

## Folder classification

| Folder | Contents and use |
| --- | --- |
| `code/` | Newly prepared Python consistency checker/exporter and base-R final-network plotting helper |
| `data/processed/` | Supplied 59 consensus targets and final 11-gene recurrence statistics |
| `results/miRNA/` | Supplied discovery and validation summaries |
| `results/network/` | Cytoscape node and edge CSVs, excluded genes, final axis, network degree tables |
| `figures/` | Original supplied network PNG; R helper can create a new schematic PDF |
| `source/original_scripts/` | Selected original R files preserved unchanged; inspect relevant blocks before execution |
| `source/original_tables/` | Original workbook exports and earlier 14-gene axis retained for provenance |
| `documentation/` | Provenance, file classification, checksums, validation report and unresolved release requirements |

## Run the processed-data checks

Python 3 with no third-party packages is sufficient:

```bash
python3 code/validate_and_export_network.py
```

This reads the supplied CSVs, checks 59 consensus genes, 11 final genes, 14 network nodes and 13 edges, verifies miRNA assignments and agreement of previous statistics, and writes:

- `results/network/final_11_gene_ceRNA_axis.csv`
- `results/network/final_ceRNA_node_degree.csv`
- `documentation/validation_report.json`

These outputs are already included. This command does not recompute expression differences or establish database provenance.

To generate a new schematic of the supplied final network using base R:

```bash
Rscript code/plot_final_network.R
```

The plotting helper was reviewed statically but was not executed here. Its output is a reconstructed schematic, not a claim to reproduce the exact original manuscript image.

## Data interpretation

The final network comprises hsa_circRNA_000554, hsa-miR-182-5p, hsa-miR-590-5p and 11 mRNAs. PDCD4 belongs to the miR-590-5p branch; the other ten genes belong to the miR-182-5p branch. CBX4, STAT3 and THBS1 are excluded from the downregulated network.

The recurrence comparison is Recurrence minus No recurrence. Selection uses nominal `P.Value < 0.05` and `logFC < 0`. All 11 genes have `adj.P.Val > 0.05`; they are exploratory nominally significant candidates, not FDR-significant genes. `adj.P.Val` and probe IDs are preserved in the final table.

`consensus_59_target_network_degree.csv` describes the two-miRNA/59-target network; it is not a STRING PPI degree table. `final_ceRNA_node_degree.csv` describes the full circRNA/miRNA/11-mRNA network; it is also not PPI topology.

The GSE154255 table supplied has 24 rows. It is a selected differential-expression export, not the entire measured expression matrix or all tested features. GSE38167 validation maps discovery hsa-miR-182-5p to validation hsa-miR-182; the corresponding arm annotation should be documented in the final methods.

## Original R sources

The strongest final-network source is `12_ceRNA_Network_Final_11_GSE43502_Sankey_Final_11_UPDATED.R`. It includes the whole earlier final-network source and an additional Sankey section. It also includes backup actions, local Windows paths and dependencies on external RDS/CSV/workspace files. Do not execute it unchanged as an automated pipeline.

`05_GSE101124_circRNA_validation.R` is retained for its target-database parsing and intersection blocks. Its GSE101124 validation section is exploratory; the manuscript uses GSE182471 for final circRNA validation.

The recurrence source contains an unfinished expression and several repeated troubleshooting sections. The circRNA validation source derives an earlier probe identifier from a different dataset while later labeling ASCRP3009378. Correct and verify extraction before asserting reproducibility. Older 14-gene network and survival scripts are excluded from active code.

Original workbooks have inconsistent declared worksheet dimensions and a missing drawing reference. Their useful tables were read using actual worksheet rows; node and edge CSVs are the preferred direct imports. Original workbook bytes are preserved without modification.

## Public datasets

GEO identifiers reported in the manuscript are GSE154255, GSE38167, GSE101123, GSE182471 and GSE43502. Sample expression and phenotype inputs for GSE43502 were subsequently recovered from the supplied RDS objects and are now included. Inputs for the other datasets remain incomplete. Original scripts contain GEO download operations, but those operations were not run or verified during release preparation.

## Citation and release metadata

No DOI, public repository URL, license or software version history is invented in this package. Add actual repository/release identifiers after publication. Choose a license before releasing code; respect third-party database terms for their exported data.

After publication of this scoped package, appropriate wording is:

“Selected analysis scripts and processed data for the candidate ceRNA network are available at [actual repository or archive DOI]. Public transcriptomic datasets are identified by their GEO accession numbers in Table 1.”

Do not use “all analyses are fully reproducible” until the missing material has been recovered and the complete pipeline verified.

## Added STRING PPI and ROC material

The original uploaded STRING query contains 14 genes and 17 interactions, not the final 11-gene set. It is preserved in `source/original_ppi_14_genes/`. The included script filters that export to the final 11 genes and combined_score >= 0.4:

```bash
python3 code/filter_and_validate_ppi.py
```

This produces nine undirected edges among eight connected genes. BCL2 has degree 5; FOXO3 and PDCD4 have degree 3 each. CYLD, LSM14A and PPM1L are isolated. These counts match the core membership and hub described in the manuscript. This is an induced subgraph of the supplied prior query, not an independently queried final STRING export; STRING version/access date were not supplied. Original workbooks named FINAL have BCL2 degree 7 because they describe the earlier 14-gene analysis; use the derived final PPI CSVs for the current network. All original edge scores and evidence channels are preserved in the filtered CSV.

The uploaded GSE182471 ROC summary is in `results/circRNA/`. It reports AUC 0.96 and confidence interval 0.849127694052026 to 1, with higher expression in normal than tumor. It was copied unchanged. Confidence level/method, sample expression/group data and ROC implementation were not supplied, so the AUC and interval were not recomputed.

## Added TCGA subtype boxplot summaries

Eleven supplied UALCAN-style CSV exports are preserved in `source/UALCAN_exports/`. They contain Normal, Luminal, HER2Pos and seven TNBC subgroup five-number summaries. They are not individual patient expression measurements. Gene labels were inferred from supplied filenames. The source CSVs do not state units, platform version, query access date or statistical comparison P values. The manuscript calls the values TPM; confirm that interpretation against the original query.

```bash
python3 code/combine_ualcan_exports.py
```

This produces `results/UALCAN/subtype_expression_boxplot_summary.csv` with 110 rows and checks group counts and low <= Q1 <= median <= Q3 <= high. Present genes: BCL2, CADM1, CASP2, CYLD, FOXO3, HOXA9, KDM5A, LSM14A, PDCD4, PPM1L and RECK. All 11 final genes have subgroup boxplot summaries. The two uploaded RECK copies are byte-identical; one is included. The combined file preserves low/high as exported, without assuming they are observed minima/maxima rather than whiskers. No significance tests are inferred from these summaries.

Manuscript reconciliation: supplied PDCD4 median expression is 55.97 in TNBC-M, 101.92 in TNBC-LAR and 76.72 in TNBC-BL1. The statement of relatively higher PDCD4 in M needs review against these actual exports and a specified reference group. Luminal is one combined category in these CSVs, not a separate Luminal A category.

## Supplied PPM1L survival plot

`source/earlier_survival/` contains the original uploaded PPM1L Kaplan-Meier image and transcribed statistics (HR 1.18, displayed CI 0.65–2.15, log-rank P 0.59). The filename labels RFS; the image does not identify endpoint or cohort filters. This older-source figure is not substituted for the manuscript overall-survival analysis and is separate from the now-supplied UALCAN PPM1L expression export.

## Added final GO and KEGG exports

`results/enrichment/` preserves the supplied final 11-gene Entrez mapping, GO CSV, KEGG CSV, workbook and no-significant-KEGG explanation. The GO export has 144 terms with adjusted P < 0.05. The KEGG file intentionally contains column headers and zero result rows; it is not a missing or accidentally blank file, according to the supplied analysis note.

```bash
python3 code/validate_enrichment_tables.py
```

The checker verifies all 11 mapped genes, term gene membership, Count/GeneRatio agreement, valid P/q values and export counts. These are validation checks on supplied results, not rerun enrichment. The saved effective annotation universes and settings were subsequently recovered from enrichment RDS objects. Custom background choice, annotation versions and executed R session remain unverified.

## Final supplied-material release preparation

Repository/release descriptions and a scoped manuscript availability statement are provided in `documentation/PUBLICATION_TEXT.md`. The newly uploaded `.RData(3)` is byte-identical to the earlier workspace and could not be decoded by the available reader; no R runtime was available. This is a reader failure, not proof of corruption. No unverified workspace contents have been added. See `documentation/R_workspace_assessment.json`.

## Supplied R command history

The supplied 512-line history is archived unchanged in `source/original_history/`. It contains GSE154255 plotting and export commands, but no ROC or other missing analysis commands. It is incomplete and depends on pre-existing R objects; it is documentation, not a runnable pipeline. Its accompanying README lists candidate output filenames for later recovery.

## Original backup description

The supplied original backup README is preserved in `source/original_backup_documentation/`. Its final gene sets and PPI hub agree with the supplied tables. It names the local backup folder `99_FINAL_PROJECT_BACKUP_TNBC_ceRNA_Recurrence` and describes additional RDS objects and a manifest; those objects are not established as supplied by the description alone. The second 512-line history contains exploratory GSE45498/GSE41970 commands, not the missing ROC implementation, and is excluded from active manuscript code.

## Recovered original RDS objects and recurrence data

All 13 supplied RDS objects are preserved unchanged in `source/original_rds/`. `data/recurrence/` provides accessible exports, including a gzip CSV expression matrix (54,675 probes × 25 samples), full probe-level results, sample metadata and sample group mapping. These are supplied analysis-stage expression values, not a claim of raw array intensities. Sample order agrees with GSM identifiers and groups: 16 Recurrence and 9 No_recurrence. All probe logFC values agree with the supplied group mean differences (maximum absolute discrepancy approximately 7.1e-15). Final 11-gene logFC/P/adjusted-P, final network edges/nodes and GO table agree with previous exports. Limma t/P/FDR statistics were not independently recomputed.

The `enrichResult` objects were decoded into their underlying slots, without executing R/Bioconductor methods. Exported settings are BH, pvalueCutoff 0.05 and qvalueCutoff 0.2; effective annotation universes contain 18,860 GO and 9,399 KEGG Entrez IDs. Saved term gene sets are included. The KEGG object stores 62 computed rows, none with adjusted P < 0.05; the final filtered dataframe remains empty. The stored 62-row table is not a significant-pathway export. Annotation/package versions and whether a custom background was chosen are not established by these objects.

To regenerate the recovered exports and their validation report, install Python dependencies `rdata`, `pandas`, `numpy` and `xarray`, then run:

```bash
python3 code/export_and_validate_rds.py
```

This exporter was executed successfully during preparation. It validates/exports stored objects and does not rerun differential expression or enrichment. Final circRNA ROC input/code and final OS results remain absent.
