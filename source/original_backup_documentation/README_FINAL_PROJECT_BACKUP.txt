TNBC ceRNA recurrence project - FINAL BACKUP
============================================

Main analysis:
- hsa_circRNA_000554-centered ceRNA network
- miRNAs: hsa-miR-182-5p and hsa-miR-590-5p
- Validation dataset: GSE43502
- Comparison: Recurrence vs No_recurrence

Final gene sets:
- 59 common miRNA target genes from TargetScan, miRDB, and miRTarBase
- 11 validated downregulated recurrence-associated ceRNA target genes
- 8 STRING-supported PPI core genes
- Hub gene in PPI: BCL2

Final 11 downregulated genes:
FOXO3, RECK, BCL2, PPM1L, KDM5A, CADM1, CYLD, HOXA9, CASP2, LSM14A, PDCD4

PPI core genes:
BCL2, FOXO3, PDCD4, RECK, CADM1, CASP2, HOXA9, KDM5A

Excluded from final downregulated ceRNA network due to upregulation:
CBX4, STAT3, THBS1

Genes not included in PPI core due to lack of STRING interactions above score 0.4:
PPM1L, CYLD, LSM14A

GO/KEGG:
- GO enrichment was performed for the final 11 genes.
- KEGG enrichment did not identify significant pathways for the final 11 genes.

Backup contains:
- R workspace
- RDS objects
- ceRNA network files
- PPI files
- GO/KEGG files
- figures
- Excel/CSV tables
- project file manifest

Backup folder:
99_FINAL_PROJECT_BACKUP_TNBC_ceRNA_Recurrence 

Created on:
2026-06-23 02:33:19.153009 
