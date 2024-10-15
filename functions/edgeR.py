# edgeR.py

#under construction
#import numpy as np
#import pandas as pd
#from scipy import stats
#from statsmodels.formula.api import glm
#import statsmodels.api as sm


#EdgeR recreated in python

def calc_norm_factors(counts):
    """
    Calculate normalization factors using the TMM method.
    """
    # Reference sample (use the sample with median library size)
    lib_sizes = counts.sum()
    ref_column = counts.columns[np.argmin(abs(lib_sizes - np.median(lib_sizes)))]
    ref_counts = counts[ref_column]

    norm_factors = {}
    for sample in counts.columns:
        if sample == ref_column:
            norm_factors[sample] = 1.0
            continue
        test_counts = counts[sample]
        log_ratios = np.log2((test_counts + 1) / (ref_counts + 1))
        abs_log_ratios = np.abs(log_ratios)
        sum_counts = test_counts + ref_counts

        # Exclude genes with low counts or extreme log ratios
        valid = (sum_counts > 0) & (abs_log_ratios < np.log2(10))
        trimmed = log_ratios[valid]
        # Trim the log ratios (e.g., remove top and bottom 30%)
        trim_fraction = 0.3
        lower = trimmed.quantile(trim_fraction)
        upper = trimmed.quantile(1 - trim_fraction)
        trimmed = trimmed[(trimmed >= lower) & (trimmed <= upper)]
        # Compute the mean of trimmed log ratios
        norm_factors[sample] = 2 ** (trimmed.mean())
    # Normalize so that the product of norm factors equals 1
    scaling_factor = np.prod(list(norm_factors.values())) ** (1 / len(norm_factors))
    norm_factors = {k: v / scaling_factor for k, v in norm_factors.items()}
    return pd.Series(norm_factors)

def estimate_common_dispersion(counts, groups, norm_factors):
    """
    Estimate common dispersion across all genes.
    """
    lib_sizes = counts.sum()
    normed_lib_sizes = lib_sizes * norm_factors
    normed_counts = counts.div(norm_factors, axis=1)

    group_labels = np.unique(groups)
    if len(group_labels) != 2:
        raise ValueError("Common dispersion estimation is implemented for two groups only.")

    group1 = counts.columns[groups == group_labels[0]]
    group2 = counts.columns[groups == group_labels[1]]

    dispersions = []
    for gene in counts.index:
        y1 = normed_counts.loc[gene, group1]
        y2 = normed_counts.loc[gene, group2]
        mean1 = y1.mean()
        mean2 = y2.mean()
        var1 = y1.var()
        var2 = y2.var()
        # Pooled variance estimate
        pooled_var = ((len(y1) - 1) * var1 + (len(y2) - 1) * var2) / (len(y1) + len(y2) - 2)
        if mean1 > 0 and mean2 > 0:
            dispersion = (pooled_var - (mean1 + mean2) / 2) / ((mean1 + mean2) / 2) ** 2
            dispersions.append(dispersion)
    common_dispersion = np.median(dispersions)
    return max(common_dispersion, 1e-8)  # Avoid negative or zero dispersion

def estimate_tagwise_dispersion(counts, groups, norm_factors, common_dispersion):
    """
    Estimate tagwise (gene-wise) dispersion.
    """
    # For simplicity, we'll shrink gene-wise dispersions towards the common dispersion
    tagwise_dispersion = {}
    for gene in counts.index:
        y = counts.loc[gene]
        lib_sizes = counts.sum()
        normed_lib_sizes = lib_sizes * norm_factors
        mu = y.sum() / normed_lib_sizes.sum() * normed_lib_sizes
        variance = mu + mu ** 2 * common_dispersion
        residual = (y - mu) ** 2 - variance
        dispersion = (residual.sum()) / (mu ** 2).sum()
        # Shrink towards common dispersion
        dispersion = 0.1 * dispersion + 0.9 * common_dispersion
        tagwise_dispersion[gene] = max(dispersion, 1e-8)
    return pd.Series(tagwise_dispersion)

def exact_test(counts, groups, norm_factors, dispersions):
    """
    Perform exact tests for differential expression between two groups.
    """
    group_labels = np.unique(groups)
    if len(group_labels) != 2:
        raise ValueError("Exact test is implemented for two groups only.")

    group1 = counts.columns[groups == group_labels[0]]
    group2 = counts.columns[groups == group_labels[1]]

    results = []
    for gene in counts.index:
        y1 = counts.loc[gene, group1]
        y2 = counts.loc[gene, group2]
        n1 = len(y1)
        n2 = len(y2)
        dispersion = dispersions[gene]
        size_factors1 = norm_factors[group1]
        size_factors2 = norm_factors[group2]
        lib_sizes1 = counts[group1].sum() * size_factors1
        lib_sizes2 = counts[group2].sum() * size_factors2
        # Compute the exact test statistic (for simplicity, we'll use a normal approximation)
        mu1 = y1.mean()
        mu2 = y2.mean()
        var1 = mu1 + mu1 ** 2 * dispersion
        var2 = mu2 + mu2 ** 2 * dispersion
        se = np.sqrt(var1 / n1 + var2 / n2)
        if se == 0:
            p_value = 1.0
            log_fc = 0.0
        else:
            z = (mu1 - mu2) / se
            p_value = 2 * (1 - stats.norm.cdf(abs(z)))
            log_fc = np.log2((mu2 + 1e-8) / (mu1 + 1e-8))
        results.append({
            'gene': gene,
            'logFC': log_fc,
            'p_value': p_value
        })
    results_df = pd.DataFrame(results).set_index('gene')
    return results_df

def glm_lrt(counts, design, norm_factors, dispersions, coef_name):
    """
    Fit a negative binomial GLM and perform likelihood ratio tests.
    """
    results = []
    for gene in counts.index:
        counts_gene = counts.loc[gene]
        # Offset by log(norm_factors)
        offset = np.log(norm_factors)
        # Create dataframe for modeling
        data = design.copy()
        data['counts'] = counts_gene
        data['offset'] = offset
        # Fit full model
        nb_family = sm.families.NegativeBinomial(alpha=dispersions[gene])
        formula_full = 'counts ~ ' + ' + '.join(design.columns[1:])  # Adjust as needed
        model_full = glm(formula_full, data=data, family=nb_family, offset=data['offset']).fit()
        # Fit reduced model (without the coefficient of interest)
        reduced_terms = [term for term in design.columns[1:] if term != coef_name]
        formula_reduced = 'counts ~ ' + ' + '.join(reduced_terms) if reduced_terms else 'counts ~ 1'
        model_reduced = glm(formula_reduced, data=data, family=nb_family, offset=data['offset']).fit()
        # Likelihood ratio test
        lr_stat = 2 * (model_full.llf - model_reduced.llf)
        p_value = stats.chi2.sf(lr_stat, df=1)
        log_fc = model_full.params.get(coef_name, np.nan) / np.log(2)
        results.append({
            'gene': gene,
            'logFC': log_fc,
            'p_value': p_value
        })
    results_df = pd.DataFrame(results).set_index('gene')
    return results_df
