"""
## dependancy
You should install pandas firstly:

$ pip install pandas

## About FPKM
According to GDC's definition:
(https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-expression-normalization)

FPKM = (RCg * 10**9) / (RCpc * L)
FPKM-UQ = (RCg * 10**9) / (RCg75 * L)

---------------
RCg: Number of reads mapped to the gene
RCpc: Number of reads mapped to all protein-coding genes
RCg75: The 75th percentile read count value for genes in the sample
L: Length of the gene in base pairs; Calculated as the sum of all exons in a gene
---------------

"""

import re
import sys
import argparse

import pandas as pd


parser = argparse.ArgumentParser(
    description="script for merged htseq-count table(tab delimated) normalization.")
parser.add_argument("table")
parser.add_argument("gtf")
parser.add_argument("--method", "-m",
    default="FPKM",
    help="The method of normalization.")


def get_gene_id(line):
    match = re.match('.*gene_id \"(.*?)\";', line)
    return match.groups()[0]

if __name__ == "__main__":
    args = parser.parse_args()

    # parse gtf
    protein_coding_genes = set([])
    gene_length = {}

    sys.stderr.write("parsing gtf...\n")
    sys.stderr.flush()

    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#"): # comment line
                continue
            items = line.split(sep="\t")
            name = get_gene_id(line)
            start, end = int(items[3]), int(items[4])
            tp = items[2]
            if tp == 'gene': # record is gene
                if "protein_coding" in line:
                    protein_coding_genes.add(name)
                gene_length[name] = abs(start - end)


    sys.stderr.write("normalization...\n")
    sys.stderr.flush()

    counts = pd.read_csv(args.table, sep='\t', index_col=0)
    pc_counts = pd.DataFrame( # subset of counts for protein coding genes
        index=pd.Series(list(protein_coding_genes))).join(counts)

    RCpc = pc_counts.sum()
    RCg75 = pc_counts.quantile(0.75)

    genes, lengthes = zip(*gene_length.items())
    L = pd.DataFrame({col:lengthes for col in list(counts)}, index=genes) # gene length table

    RCg = pd.DataFrame(index=L.index).join(counts) # subset of counts filter non-gene rows

    if args.method == "FPKM":
        res = (RCg * 10**9) / (RCpc * L)
    elif args.method == "FPKM-UQ":
        res = (RCg75 * 10**9) / (RCg75 * L)
    else:
        raise NotImplementedError(args.method + " method have not implemented.")

    res.index.name = 'gene'

    print(res.to_csv(sep="\t"))
