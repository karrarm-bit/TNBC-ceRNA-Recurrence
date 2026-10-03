# Supplied R history

GSE154255_supplied_Rhistory.txt preserves all 512 supplied history lines unchanged. It is an incomplete interactive history, starting within a ggplot expression, and is not a standalone executable R script. It records PCA, volcano, boxplot and heatmap commands and CSV/figure export commands associated with GSE154255. No ROC, GSE182471, GSE43502, GO/KEGG enrichment or survival commands were found. Existing objects such as expr_tnbc, group and deg are referenced without their construction in this history. A history contains commands, not their outputs or evidence of successful execution.

The export commands reference the following potentially recoverable local results:
- results/01_all_differential_expression_results.csv
- results/02_upregulated_miRNAs_Pvalue_only.csv
- results/03_downregulated_miRNAs_Pvalue_only.csv
- results/04_selected_2_upregulated_2_downregulated_miRNAs.csv
- results/05_expression_values_hsa_miR_590_5p.csv (generated if the matched row name was exactly hsa-miR-590-5p)
- results/05_expression_values_hsa_miR_182_5p.csv (generated if the matched row name was exactly hsa-miR-182-5p)
- results/06_top20_miRNAs_for_heatmap.csv

These paths are recorded commands; existence and file contents were not verified. This history does not recover ROC implementation or missing sample data.
