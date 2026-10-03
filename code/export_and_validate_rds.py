"""Export supplied RDS using rdata, pandas, numpy and xarray; does not rerun limma."""
from pathlib import Path
import json,shutil
import numpy as np
import pandas as pd
import rdata
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'source/original_rds'
OUT=ROOT/'data/recurrence';OUT.mkdir(parents=True,exist_ok=True)
objs={p.stem:rdata.read_rds(p) for p in sorted(SRC.glob('*.rds'))}
x=objs['expr_prognosis'];ph=objs['pheno_prognosis'];groups=np.asarray(objs['group_mrna']).astype(str)
samples=x.coords[x.dims[1]].values.astype(str);probes=x.coords[x.dims[0]].values.astype(str)
assert x.shape==(54675,25) and np.isfinite(x.values).all()
assert list(samples)==list(ph['GSM'].astype(str))
assert list(groups)==list(ph['Prognosis'].astype(str))
assert dict(pd.Series(groups).value_counts())=={'Recurrence':16,'No_recurrence':9}
ph.to_csv(OUT/'GSE43502_sample_metadata.csv',index=False)
pd.DataFrame({'GSM':samples,'Group':groups}).to_csv(OUT/'GSE43502_sample_groups.csv',index=False)
pd.DataFrame(x.values,index=pd.Index(probes,name='ProbeID'),columns=samples).to_csv(OUT/'GSE43502_supplied_expression_matrix.csv.gz',compression='gzip')
mr=objs['mrna_results'];assert len(mr)==54675 and mr['ProbeID'].is_unique
mr.to_csv(OUT/'GSE43502_full_mrna_results.csv.gz',index=False,compression='gzip')
old=pd.read_csv(ROOT/'data/processed/GSE43502_final_11_downregulated_targets.csv')
print('final_table_columns',list(old.columns))
probe_col='ProbeID';matches=[]
for _,row in old.iterrows():
 hit=mr.loc[mr[probe_col].astype(str)==str(row[probe_col])];assert len(hit)==1
 for c in ['logFC','P.Value','adj.P.Val']:
  assert np.isclose(float(hit.iloc[0][c]),float(row[c]),rtol=1e-8,atol=1e-12),(row[probe_col],c)
 matches.append(str(row[probe_col]))
idx={p:i for i,p in enumerate(probes)}
means=x.values[:,groups=='Recurrence'].mean(axis=1)-x.values[:,groups=='No_recurrence'].mean(axis=1)
diffs=[abs(means[idx[str(row.ProbeID)]]-float(row.logFC)) for row in mr.itertuples()]
assert max(diffs)<1e-7
final=set(map(str,objs['final_downregulated_genes']));assert len(final)==11 and final==set(old['Gene'])
for name in ['final_nodes','final_edges','final_11_network_genes','gene_df_11','ego_11_df','ekegg_11_df']:
 objs[name].to_csv(OUT/(name+'_from_RDS.csv'),index=False)
assert set(objs['final_nodes']['name'])==set(pd.read_csv(ROOT/'results/network/ceRNA_nodes.csv')['name'])
e=objs['final_edges'];prior=pd.read_csv(ROOT/'results/network/ceRNA_edges.csv');assert set(zip(e.source,e.target))==set(zip(prior.source,prior.target))
assert set(objs['final_11_network_genes']['Gene'])==final
pd.testing.assert_frame_equal(objs['ego_11_df'].reset_index(drop=True),pd.read_csv(ROOT/'results/enrichment/GO_enrichment_final_11_genes.csv'),check_dtype=False,check_exact=False,rtol=1e-8)
ENR=ROOT/'results/enrichment';settings={}
for n in ['ego_11','ekegg_11']:
 o=objs[n];meta={k:np.asarray(getattr(o,k)).tolist() for k in ['pvalueCutoff','pAdjustMethod','qvalueCutoff','organism','ontology','gene','keytype','readable']}
 meta['effective_universe_size']=len(o.universe);meta['stored_result_rows']=len(o.result);meta['adjusted_P_lt_0_05']=int((o.result['p.adjust']<.05).sum());settings[n]=meta
 pd.DataFrame({'EntrezID':np.asarray(o.universe)}).to_csv(ENR/(n+'_effective_universe.csv'),index=False)
 o.result.to_csv(ENR/(n+'_stored_result_table.csv'),index=False)
 with (ENR/(n+'_stored_gene_sets.json')).open('w') as f:json.dump({str(k):np.asarray(v).tolist() for k,v in o.geneSets.items()},f)
(ENR/'recovered_enrichment_settings.json').write_text(json.dumps(settings,indent=2)+'\n')
report={'status':'passed','RDS_files_decoded':len(objs),'expression_shape':list(x.shape),'sample_order_matches_metadata':True,'group_counts':{'Recurrence':16,'No_recurrence':9},'all_probe_logFC_matches_group_mean_difference':True,'max_logFC_difference':max(diffs),'final_11_statistics_match_prior_export':True,'final_network_membership_and_edges_match':True,'GO_dataframe_matches_prior_export':True,'GO_effective_universe':len(objs['ego_11'].universe),'KEGG_effective_universe':len(objs['ekegg_11'].universe),'KEGG_stored_rows':len(objs['ekegg_11'].result),'KEGG_adjusted_significant_rows':0,'scope':'RDS decoding, export and consistency checks; limma p-values and enrichment were not independently recomputed; effective annotation universes do not establish a user-specified background or database versions.'}
(ROOT/'documentation/rds_validation_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
