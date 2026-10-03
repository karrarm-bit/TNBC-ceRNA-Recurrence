"""Normalize supplied boxplot summaries; does not create individual samples or P values."""
import csv,re,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 rows=[];genes=set();group_counts={}
 for p in sorted((ROOT/'source/UALCAN_exports').glob('expression-of-*.csv')):
  gene=re.search(r'expression-of-([^-]+)-in',p.name).group(1).upper();genes.add(gene)
  with p.open(encoding='utf-8-sig',newline='') as f: records=list(csv.DictReader(f))
  if len(records)!=10:raise ValueError('Expected 10 subgroup summary rows: '+p.name)
  for r in records:
   m=re.fullmatch(r'(.*?)<br>\(n=(\d+)\)',r['TCGA samples'])
   if not m:raise ValueError('Unrecognized group label')
   group,n=m.group(1),int(m.group(2))
   if group in group_counts and group_counts[group]!=n:raise ValueError('Different sample counts')
   group_counts[group]=n
   values=[float(r['Series 1 ('+k+')']) for k in ['low','q1','median','q3','high']]
   if not all(math.isfinite(v) for v in values) or values!=sorted(values):raise ValueError('Invalid boxplot statistics')
   rows.append(dict(zip(['Gene','Group','N','Low','Q1','Median','Q3','High','Source_file'],[gene,group,n,*values,p.name])))
 with (ROOT/'results/UALCAN/subtype_expression_boxplot_summary.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 expected={'FOXO3','PDCD4','RECK','BCL2','PPM1L','KDM5A','CADM1','CYLD','HOXA9','CASP2','LSM14A'}
 report={'status':'passed','files':len(genes),'summary_rows':len(rows),'genes_present':sorted(genes),'missing_final_genes':sorted(expected-genes),'group_counts':group_counts,'data_grain':'One gene/subgroup boxplot summary, not individual samples','P_values_available':False,'units':'Not present in CSV; manuscript states TPM, verify platform export settings.'}
 (ROOT/'documentation/ualcan_validation_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
