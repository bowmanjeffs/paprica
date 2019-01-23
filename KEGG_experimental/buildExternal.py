import os
import sys

ref_dir = 'ref_genome_database/'
model_dir='models/'
def alignDomain(domain):
    cwd=os.getcwd()+'/'
    ref_dir_domain=cwd+ref_dir+domain+'/'
    model_dir_domain=model_dir+domain+'_ssu.cm'
    align_dir=ref_dir_domain+'align/'
    os.mkdir(align_dir)
    os.system('cmalign --dna --outformat Pfam -o '+align_dir+domain+'.align.DNA.sto '+model_dir_domain+' '+ref_dir_domain+'16S_'+domain+'_kegg.fasta')
    os.system('seqmagick convert '+align_dir+domain+'.align.DNA.sto '+align_dir+domain+'.align.DNA.fasta')

def treeDomain(domain):
    cwd=os.getcwd()+'/'
    ref_dir_domain=cwd+ref_dir+domain+'/'
    model_dir_domain=model_dir+domain+'_ssu.cm'
    align_dir=ref_dir_domain+'align/'
    tree_dir=ref_dir_domain+'tree/'
    os.mkdir(tree_dir)
    os.system('raxmlHPC-PTHREADS -T 8 -m GTRGAMMA -s '+align_dir+domain+'.align.DNA.fasta -f d -p 12345 -n '+domain+'.ref.tre -w '+tree_dir)
    os.system('raxmlHPC-PTHREADS -T 2 -m GTRGAMMA -t '+tree_dir+'RAxML_bestTree.'+domain+'.ref.tre -f I -p 12345 -n '+domain+'.root.ref.tre -w '+tree_dir)
    os.system('raxmlHPC-PTHREADS -T 8 -m GTRGAMMA -f J -p 12345 -t '+tree_dir+'RAxML_rootedTree.'+domain+'.root.ref.tre -s '+align_dir+domain+'.align.DNA.fasta -n '+domain+'.conf.root.ref.tre -w '+tree_dir)

def taxtableDB(domain):
    cwd=os.getcwd()+'/'
    ref_dir_domain=cwd+ref_dir+domain+'/'
    os.system('taxit taxtable '+ref_dir+'ncbi_taxonomy.db -f '+ref_dir_domain+'tax_ids.txt -o '+ref_dir_domain+'taxa.csv')


def taxtasticPackage(domain):
    cwd=os.getcwd()+'/'
    ref_dir_domain=ref_dir+domain+'/'
    model_dir_domain=model_dir+domain+'_ssu.cm'
    align_dir=ref_dir_domain+'align/'
    tree_dir=ref_dir_domain+'tree/'
    os.system('taxit create -l 16S_rRNA -P '+ref_dir_domain+'taxit ' \
					'--aln-fasta '+align_dir+domain+'.align.DNA.fasta ' \
					'--aln-sto '+align_dir+domain+'.align.DNA.sto ' \
					'--tree-file '+tree_dir+'RAxML_fastTreeSH_Support.'+domain+'.conf.root.ref.tre ' \
					'--tree-stats '+tree_dir+'RAxML_info.'+domain+'.ref.tre ' \
					'--seq-info '+ref_dir_domain+'seq_info.csv ' \
					'--taxonomy '+ref_dir_domain+'taxa.csv' \
		)

def main():
    for domain in ["archaea","bacteria"]:
        alignDomain(domain)
        treeDomain(domain)
        taxtableDB(domain)
        taxtasticPackage(domain)

if __name__ == "__main__":
    main()
