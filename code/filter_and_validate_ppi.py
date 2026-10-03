"""Derive final 11-gene PPI from the supplied earlier 14-gene STRING export.
No network calls; Python standard library only. Scores preserved from supplied export.
"""
from pathlib import Path
import csv,json
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
def read(p,delimiter=','):
    with (ROOT/p).open(newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))
def write(p,rows,fields):
    with (ROOT/p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def main():
    original=read('source/original_ppi_14_genes/string_interactions_short.tsv','\t')
    nodes=read('source/original_ppi_14_genes/string_node_degrees.tsv','\t')
    original_degree=Counter()
    for r in original:original_degree.update([r['#node1'],r['node2']])
    if len(original)!=17 or any(int(r['node_degree'])!=original_degree[r['#node']] for r in nodes):
        raise ValueError('Original STRING degree/edge exports disagree')
    genes={r['Gene'] for r in read('data/processed/GSE43502_final_11_downregulated_targets.csv')}
    final=[r for r in original if r['#node1'] in genes and r['node2'] in genes and float(r['combined_score'])>=0.4]
    pairs={frozenset([r['#node1'],r['node2']]) for r in final}
    if len(pairs)!=len(final):raise ValueError('Repeated undirected edges')
    degree=Counter()
    for r in final:degree.update([r['#node1'],r['node2']])
    expected={'BCL2','FOXO3','PDCD4','RECK','CADM1','CASP2','HOXA9','KDM5A'}
    if set(degree)!=expected or len(final)!=9 or degree['BCL2']!=5:
        raise ValueError('Filtered PPI does not match expected manuscript core')
    write('results/PPI/final_11_gene_STRING_edges.csv',final,list(original[0]))
    rows=[{'Gene':g,'Degree':degree[g],'Status':'Connected' if degree[g] else 'Isolated'} for g in sorted(genes,key=lambda g:(-degree[g],g))]
    write('results/PPI/final_11_gene_PPI_degree.csv',rows,['Gene','Degree','Status'])
    write('results/PPI/final_11_gene_PPI_isolated.csv',[r for r in rows if not r['Degree']],['Gene','Degree','Status'])
    report={'status':'passed','original_query_genes':14,'original_edges':17,'final_input_genes':11,'final_connected_genes':8,'final_edges':9,'threshold':0.4,'hub':'BCL2','hub_degree':5,'isolated':['CYLD','LSM14A','PPM1L'],'provenance':'Induced 11-gene subgraph filtered from supplied 14-gene STRING export; not a new STRING query.'}
    (ROOT/'documentation/ppi_validation_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
