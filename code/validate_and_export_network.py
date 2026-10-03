"""Validate supplied tables and export the final ceRNA axis. Python standard library only.
Run from any directory: python3 code/validate_and_export_network.py
This validates supplied processed results; it does not rerun GEO differential expression.
"""
from pathlib import Path
import csv, json, math
from collections import Counter
ROOT = Path(__file__).resolve().parents[1]
def read(path):
    with (ROOT/path).open(newline='',encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))
def write(path,rows,fields):
    with (ROOT/path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def require(condition,message):
    if not condition: raise ValueError(message)
def main():
    targets=read('data/processed/consensus_miRNA_targets_59.csv')
    genes=read('data/processed/GSE43502_final_11_downregulated_targets.csv')
    nodes=read('results/network/ceRNA_nodes.csv')
    edges=read('results/network/ceRNA_edges.csv')
    old=read('source/original_tables/earlier_14_gene_ceRNA_axis.csv')
    degree=read('results/network/consensus_59_target_network_degree.csv')
    expected={'FOXO3','PDCD4','RECK','BCL2','PPM1L','KDM5A','CADM1','CYLD','HOXA9','CASP2','LSM14A'}
    require(len(targets)==59 and len({r['Gene'] for r in targets})==59,'Expected 59 unique consensus genes')
    require(Counter(r['miRNA'] for r in targets)=={'hsa-miR-182-5p':50,'hsa-miR-590-5p':9},'Consensus branch counts differ')
    require(len(genes)==11 and {r['Gene'] for r in genes}==expected,'Final genes differ')
    require(all(float(r['logFC'])<0 and float(r['P.Value'])<0.05 for r in genes),'Final selection rule fails')
    require(all(float(r['adj.P.Val'])>=0.05 for r in genes),'Expected exploratory nominal-only selection')
    target_pairs={(r['miRNA'],r['Gene']) for r in targets}
    require(all((r['miRNA'],r['Gene']) in target_pairs for r in genes),'Final target absent from consensus')
    old_by_gene={r['Gene']:r for r in old}
    require(len(old)==14,'Earlier axis must contain 14 rows')
    for r in genes:
        previous=old_by_gene[r['Gene']]
        require(r['miRNA']==previous['miRNA'],'miRNA assignment changed')
        for col in ['logFC','P.Value']:
            require(math.isclose(float(r[col]),float(previous[col]),rel_tol=1e-7,abs_tol=1e-9),'Earlier result differs: '+r['Gene']+' '+col)
    require(len(nodes)==14 and len({r['name'] for r in nodes})==14,'Expected 14 unique nodes')
    require(Counter(r['type'] for r in nodes)=={'circRNA':1,'miRNA':2,'mRNA':11},'Node classes differ')
    node_names={r['name'] for r in nodes}
    expected_edges={('hsa_circRNA_000554','hsa-miR-182-5p'),('hsa_circRNA_000554','hsa-miR-590-5p')}|{(r['miRNA'],r['Gene']) for r in genes}
    require(len(edges)==13 and {(r['source'],r['target']) for r in edges}==expected_edges,'Final edge list differs')
    require(all(r['source'] in node_names and r['target'] in node_names for r in edges),'Edge endpoint absent')
    consensus_degree=Counter()
    for r in targets:consensus_degree.update([r['miRNA'],r['Gene']])
    require(len(degree)==61 and all(int(r['Degree'])==consensus_degree[r['Node']] for r in degree),'Consensus degrees differ')
    final_degree=Counter()
    for r in edges:final_degree.update([r['source'],r['target']])
    write('results/network/final_ceRNA_node_degree.csv',[{'Node':r['name'],'Type':r['type'],'Degree':final_degree[r['name']]} for r in nodes],['Node','Type','Degree'])
    axis=[{'circRNA':'hsa_circRNA_000554',**r} for r in genes]
    write('results/network/final_11_gene_ceRNA_axis.csv',axis,['circRNA']+list(genes[0]))
    report={'status':'passed','consensus_genes':59,'final_mRNAs':11,'final_nodes':14,'final_edges':13,'nominal_P_lt_0_05':11,'FDR_lt_0_05':0,'earlier_14_gene_values_match':True,'scope':'Consistency of supplied processed tables; no upstream GEO or web-platform analyses rerun.'}
    (ROOT/'documentation/validation_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
