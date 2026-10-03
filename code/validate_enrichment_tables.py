"""Check supplied enrichment tables and gene mappings; no enrichment is recomputed."""
from pathlib import Path
import csv,json,math
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
def read(p):
 with (ROOT/p).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def main():
 go=read('results/enrichment/GO_enrichment_final_11_genes.csv')
 kegg=read('results/enrichment/KEGG_enrichment_final_11_genes.csv')
 mapping=read('results/enrichment/Final_11_genes_EntrezIDs.csv')
 genes={r['Gene'] for r in read('data/processed/GSE43502_final_11_downregulated_targets.csv')}
 if len(mapping)!=11 or {r['SYMBOL'] for r in mapping}!=genes or len({r['ENTREZID'] for r in mapping})!=11:raise ValueError('Gene mapping mismatch')
 if len(go)!=144 or kegg:raise ValueError('Unexpected enrichment export counts')
 for r in go:
  members=r['geneID'].split('/');count=int(r['Count']);num,den=map(int,r['GeneRatio'].split('/'))
  if not set(members)<=genes or count!=len(members) or num!=count or den!=11:raise ValueError('GO gene count mismatch')
  if not all(math.isfinite(float(r[c])) and 0<=float(r[c])<=1 for c in ['pvalue','p.adjust','qvalue']):raise ValueError('Invalid P or q value')
  if float(r['p.adjust'])>=.05:raise ValueError('GO export includes adjusted P >= .05')
 report={'status':'passed','mapped_genes':11,'GO_exported_terms':144,'GO_terms_adjusted_P_lt_0_05':144,'GO_ontology_counts':dict(Counter(r['ONTOLOGY'] for r in go)),'GO_adjusted_P_range':[min(float(r['p.adjust']) for r in go),max(float(r['p.adjust']) for r in go)],'KEGG_exported_terms':0,'scope':'Validation of supplied exported results, not rerun enrichment; full tested universe and database versions not provided.'}
 (ROOT/'documentation/enrichment_validation_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
