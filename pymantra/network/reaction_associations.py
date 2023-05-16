from typing import Dict, Tuple, Union
import numpy as np
import pandas as pd
from scipy.stats import norm

from pymantra.network.enrichment.spearman import spearmans_correlation


def associate_multiomics_ld(
    residuals: Union[Dict[str, pd.Series], pd.DataFrame],
    groups: pd.Series, multiomics_data: pd.DataFrame,
    comparison: Tuple[str, str] = None,
    pval_thresh: float = 1.,
    use_pvals: bool = False
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """Associate multi-omics data with reaction changes

    Associate changes in microbiome or transcriptome compositions to changes
    in reaction activity using Spearman's correlation. More precisely
    multi-omics values are correlated with the explained variance in the linear
    reaction models (see :py:func:`pymantra.network.per_sample_ld_estimation`
    for details)

    Parameters
    ----------
    residuals : Dict[str, pd.Series] | pd.DataFrame
        (Normalized) residuals from the reaction models per reaction.
        Either a dictionary of (reaction name, pd.Series) or a
        (sample x reaction) pd.DataFrame containing the residuals
    groups : pd.Series
        Sample group annotation. The indices must match those in the values
        of `residuals`
    multiomics_data : pd.DataFrame
        Multi-omics data to associate the reactions to.
    comparison : Tuple[str, str], optional
        Groups to compare. If there are more than 2 groups in `groups`,
        this parameter is required. Else it can be used to set the order
        in which the groups are used.
    pval_thresh : float, 1.0
        p-value threshold to set correlation values to zero. With the default
        of 1 all correlation values are kept. To set all non-significant
        correlation values to 0, set this parameter to e.g. 0.05
    use_pvals : bool, False
        *Experimental*
        Whether to compute the associations against the residuals or the
        p-values of the residuals. *Currently, True is not recommended.*

    Returns
    -------
    Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]
        2-tuple where the first elements are correlations per group as a
        data frame of multi-omics x reaction and the second element are the
        correlation p-values per group in the same format as the correlations
    """
    # rows - samples, columns - reactions
    associations = {}
    pvals = {}
    # TODO: confounder correction here as well?
    if isinstance(residuals, dict):
        res_df = pd.DataFrame.from_dict(residuals)
    else:
        res_df = residuals
    # n_org = organism_data.shape[1]
    if comparison is None:
        comparison = np.unique(groups)
        if comparison.size != 2:
            raise ValueError(
                "'groups' must contain exactly 2 different groups when "
                f"'comparison' is None, but {comparison.size} unique groups "
                "found!"
            )
    for group in comparison:
        group_mask = groups == group

        n_group = group_mask.sum()
        if n_group < 3:
            raise ValueError(
                "At least three samples are required to compute multi-omics "
                f"associations, but only {n_group} found for group {group}!"
            )

        if use_pvals:
            res_pvals = 1 - (norm.sf(np.abs(res_df.loc[group_mask, :])) * 2)
            coeffs, group_pvals = spearmans_correlation(
                multiomics_data.loc[group_mask, :].values, res_pvals
            )
        else:
            coeffs, group_pvals = spearmans_correlation(
                multiomics_data.loc[group_mask, :].values,
                res_df.loc[group_mask, :].values
            )

        to_drop = np.mean(np.isnan(coeffs), axis=0) == 1

        coeffs[group_pvals > pval_thresh] = 0
        associations[group] = pd.DataFrame(
            coeffs[:, ~to_drop],
            index=multiomics_data.columns,
            columns=res_df.columns[~to_drop]
        )
        pvals[group] = pd.DataFrame(
            group_pvals[:, ~to_drop],
            index=multiomics_data.columns,
            columns=res_df.columns[~to_drop]
        )
    return associations, pvals
