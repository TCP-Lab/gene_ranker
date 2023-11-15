# Gene Ranker

This python CLI tool can rank genes based on their differential expression
between a 'case' and a 'control' status, from most up-regulated to most
down-regulated, with a variety of metrics.

This is not a differential expression analysis, as no p-values are computed since
no statistical test is used.
It is just meant to give each gene a 'rank' and a rank value in order to run
other statistical methods such as pre-ranked GSEA.

Input data should be the base-2 logarithm of 1 + read counts (i.e log2(counts + 1)).
The program un-logs the data when appropriate (e.g. running DESeq2).

Currently supported ranking methods:
- **Fold Change**: The `fold_change` method computes a simple difference of 
  average fold changes between the case and controls, without any form of
  prior normalization.
- **Normalized Cohen's D**: The `norm_cohen_d` metric normalies the data with
  the median of ratios method used by DESeq2 and implemented in pydeseq2, and
  then computes cohen's D on the resulting normalized counts.
- **DESeq2 Shrunk Log Fold Change**: Uses `DESeq2`'s LFC shrinking method to
  compute LFCs, and uses them as ranking metric.

You can use `generanker --list-methods` for a list of all the methods.

## Installation
Install Python and [`cargo`](https://doc.rust-lang.org/book/ch01-03-hello-cargo.html).
Install the prerequisite `fast-cohen` executable with:
```bash
cargo install --git https://github.com/MrHedmad/fast-cohen.git
```
Then, install the tool with:
```bash
# I suggest you do this in a virtual environment:
# python -m venv env && source env/bin/activate
python -m pip install git+https://github.com/TCP-Lab/gene_ranker.git
```
You may then use `generanker` from the command line.
Use `generanker --help` for additional usage details.

