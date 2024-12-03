from gene_ranker.methods.base import RankingMethod, norm_wrapper
from gene_ranker.methods.bws import bws_rank
from gene_ranker.methods.cohen import cohen_d_ranking
from gene_ranker.methods.deseq_shrinkage import deseq_shrinkage_ranking
from gene_ranker.methods.fold_change import fold_change_ranking
from gene_ranker.methods.signal_to_noise import signal_to_noise_ratio

RANKING_METHODS = {
    "fold_change": RankingMethod(
        name="Fold Change",
        exec=fold_change_ranking,
        parser=None,
        desc="Use a non-normalized, raw fold change metric.",
    ),
    "deseq_shrinkage": RankingMethod(
        name="DESeq2 Shrinkage",
        exec=deseq_shrinkage_ranking,
        parser=None,
        desc="Use DESeq2-shrunk fold changes. Always normalizes the input",
    ),
    "cohen_d": RankingMethod(
        name="Cohen's d",
        exec=cohen_d_ranking,
        parser=None,
        desc="Use the Cohen's d metric",
    ),
    "norm_cohen_d": RankingMethod(
        name="Normalized Cohen's d",
        exec=norm_wrapper(cohen_d_ranking),
        parser=None,
        desc="Use a DESeq2-normalized Cohen's d metric",
    ),
    "norm_fold_change": RankingMethod(
        name="Normalized Fold Change",
        exec=norm_wrapper(fold_change_ranking),
        parser=None,
        desc="Use a DESeq2-normalized fold change metric",
    ),
    "s2n_ratio": RankingMethod(
        name="Signal to noise ratio",
        exec=signal_to_noise_ratio,
        parser=None,
        desc="Use the signal to noise ratio (diff of means divided by variance)",
    ),
    "norm_s2n_ratio": RankingMethod(
        name="Normalized signal to noise ratio",
        exec=norm_wrapper(signal_to_noise_ratio),
        parser=None,
        desc="Use the signal to noise ratio metric on normalized data",
    ),
    "bws_test": RankingMethod(
        name="Baumgartner-Weiss-Schindler test statistic",
        exec=bws_rank,
        parser=None,
        desc="Use the BWS test statistic, which works well with high N samples",
    ),
    "norm_bws_test": RankingMethod(
        name="Normalized Baumgartner-Weiss-Schindler test statistic",
        exec=norm_wrapper(bws_rank),
        parser=None,
        desc="Same as BWS, but on normalized data",
    ),
}
