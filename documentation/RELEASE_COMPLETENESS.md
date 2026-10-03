# Material still required for a complete manuscript release

The final network selection table is now present and consistent. The following manuscript components are not complete in the supplied materials. Missing data have not been synthesized or replaced by assumed results.

| Component | Supplied | Still needed |
| --- | --- | --- |
| GSE154255 miRNA discovery | R source and 24-row selected results | Full tested-feature results, expression input, phenotype/group assignment and executed session/package versions |
| GSE38167 validation | R source and summary | Sample-level expression/phenotype, identifier mapping evidence and test settings |
| GSE101123 circRNA discovery | R source | Final differential table, sample-level expression/phenotype and circRNA annotation mapping |
| GSE182471 circRNA validation and ROC | R validation source and supplied ROC AUC/CI summary | Correct probe extraction evidence, sample values/groups, final validation results and actual ROC code with CI method/level; supplied ROC script is empty |
| Target intersection | 59 consensus target pairs and original parsing blocks | Actual exports from TargetScan, miRDB and miRTarBase, database versions/access dates and applied filters |
| GSE43502 recurrence | Final 11-gene statistics, all 54,675 probe results, expression matrix, 25-sample metadata/group mapping and original R source | Cleaned executed recurrence pipeline and session/package versions; p-values have not been independently rerun |
| Final ceRNA network | 14 nodes, 13 edges, final 11-gene statistics and original figure | Cytoscape session/style if exact original manuscript visualization is required |
| STRING PPI | Original 14-gene STRING edges/scores/degrees plus validated derived 11-gene subgraph, 8-gene core and final degrees | STRING version/access date/query metadata; confirm that the induced subgraph is the final network used |
| GO and KEGG | Original R blocks, final 11-gene Entrez mapping, 144 GO result rows, zero-row KEGG export, workbook and explanation | Effective annotation universes (GO 18,860; KEGG 9,399), saved gene sets and BH/p/q settings recovered from RDS. Still needed: annotation/database/package versions, executed session and clarification of custom background choice |
| Kaplan-Meier overall survival | Manuscript descriptions plus a separately archived older PPM1L plot labeled RFS in its filename, not confirmed OS | Final 11-gene platform exports and settings, HR/CI/P, probe/cut-off, group counts, cohort and access date; older forest script is incompatible |
| UALCAN subtype expression | All 11 gene/subgroup five-number summaries and a validated 110-row combined table | Statistical comparisons if claimed, original platform/unit/query metadata/access date; no individual sample data were supplied |

## Manuscript method corrections

Section 2.6 mentions adjusted P < 0.05 while section 2.7 and the actual final 11-gene table use nominal P < 0.05. Resolve this to the actual exploratory threshold. None of these 11 genes meets FDR < 0.05.

The survival methods describe R/TCGA analysis in one section and platform analysis in another, while the results describe Kaplan-Meier Plotter with an optimal cut-off. Reconcile methods with the actual platform and cohort used.

## Files deliberately excluded

- Empty placeholder scripts, including the ROC placeholder.
- Earlier 14-gene network/PPI/survival scripts whose final membership or survival values conflict with the current manuscript.
- Exploratory GSE19783, GSE121396, GSE40049, GSE45498 and GSE41970 analyses not used in the manuscript.
- The uploaded `.RData(2)` and byte-identical `.RData(3)` workspaces: it could not be decoded reliably, so its content and relevance were not verified. It is retained in the user's original upload, not distributed as validated analysis data.
- The manuscript itself: authorship/affiliations and article text do not need to be included in a code repository to reproduce calculations.

A processed-network data release can use this package with the scoped README. A complete reproducibility release for the whole manuscript needs the material listed above.

## Subtype interpretation to reconcile

PDCD4 median in the supplied TCGA-M summary is 55.97, lower than LAR 101.92, UNS 91.4, BL1 76.72 and IM 66.39. Review the manuscript wording that describes relatively higher expression in M. The exports label a combined Luminal group, not Luminal A specifically.

## RDS recovery update

All 13 supplied RDS files decoded. Recurrence sample data and full results are now supplied; their earlier absence is resolved. GO/KEGG object settings and effective universes are also supplied. See rds_validation_report.json and README for checks and remaining limits. ROC sample data/code and final manuscript OS exports were not recovered from these objects.
