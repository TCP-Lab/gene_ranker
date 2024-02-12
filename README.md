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
  average fold changes between the case and controls.
- **Cohen's d**: The `cohen_d` metric computes Cohen's d between the different
  expression values of each gene.
- **DESeq2 Shrunk Log Fold Change**: Uses `DESeq2`'s LFC shrinking method to
  compute LFCs, and uses them as ranking metric.
  This uses PyDESeq2, so the input data is always normalized in the process.
- **Signal to Noise ratio**: Compute the signal to noise ratio between the 
  control and case genes.
  This is roughly the mean divided by the variance of each gene.
- **Baumgartner-Weiss-Schindler test statistic**: Compute the
  [BWS statistic](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bws_test.html)
  for each gene. Uses a `scipy` primitive so it's much faster than using
  `bws_test`.

Most of these methods come with a normalized version, where the input is first
normalized with the ["mean of ratios" method](https://github.com/owkin/PyDESeq2/blob/39b6a373abb85991b5ac50f5f5b26a1a290d890b/pydeseq2/preprocessing.py#L8-L31)
(as implemented by DESeq2).
They are usually named as `norm_<method>`.
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

