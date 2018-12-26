# Custom scripts for RNA-Seq analyze

## Dependency

Manage dependency with `anaconda`. Download and install miniconda:

``` bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Config conda:

``` bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Install aligner, samtools, htseq:

``` bash
conda install --yes -c bioconda samtools hisat2 htseq
```

Install Python packages:

All scripts tested under Python3.6, but also try to compatible with Python2.7

```bash
conda install --yes numpy scipy pandas matplotlib seaborn scikit-learn click
conda install --yes -c bioconda pybigwig deeptools
```

Install R environment and packages:

```bash
conda install --yes -c r r-essentials
conda install --yes -c bioconda bioconductor-deseq2 bioconductor-clusterprofiler
```