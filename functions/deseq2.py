# deseq2.py

##under construction
#import numpy as np
#import pandas as pd
#from scipy import stats
#from statsmodels.formula.api import glm
#import statsmodels.api as sm

#DESeq2 recreated in python

def calculate_size_factors(counts):
    """
    Calculate size factors for normalization using the median-of-ratios method.
    """
    geometric_means = counts.apply(lambda x: stats.gmean(x + 1), axis=1)
    ratios = counts.div(geometric_means, axis=0)
    size_factors = ratios.median()
    return size_factors

def estimate_dispersion(counts, design, size_factors):
    """
    Estimate gene-wise dispersion parameters.
    """
    dispersions = {}
    for gene in counts.index:
        counts_gene = counts.loc[gene]
        # Offset by log(size factors)
        offset = np.log(size_factors)
        # Fit initial model to get mu
        data = design.copy()
        data['counts'] = counts_gene
        data['offset'] = offset
        # Fit Poisson model to get initial estimates
        formula = 'counts ~ ' + ' + '.join(design.columns[1:])  # Adjust as needed
        poisson_model = glm(formula, data=data, family=sm.families.Poisson(), offset=data['offset']).fit()
        mu = poisson_model.fittedvalues
        # Calculate dispersion (alpha)
        residuals = (counts_gene - mu)**2 - counts_gene
        alpha = max((residuals / mu**2).mean(), 1e-8)
        dispersions[gene] = alpha
    return pd.Series(dispersions)

def fit_glm_nb(counts, design, size_factors, dispersions, coef_name):
    """
    Fit a negative binomial GLM for each gene and perform hypothesis testing.
    """
    results = []
    for gene in counts.index:
        counts_gene = counts.loc[gene]
        # Log-transform size factors for offset
        offset = np.log(size_factors)
        # Create dataframe for modeling
        data = design.copy()
        data['counts'] = counts_gene
        data['offset'] = offset
        # Define the model formula
        formula = 'counts ~ ' + ' + '.join(design.columns[1:])  # Adjust as needed
        # Fit the model
        nb_family = sm.families.NegativeBinomial(alpha=dispersions[gene])
        model = glm(formula, data=data, family=nb_family, offset=data['offset']).fit()
        # Extract results for the coefficient of interest
        if coef_name in model.params:
            beta = model.params[coef_name]
            std_err = model.bse[coef_name]
            wald_stat = beta / std_err
            p_value = 2 * (1 - stats.norm.cdf(abs(wald_stat)))
            results.append({
                'gene': gene,
                'log2FoldChange': beta / np.log(2),  # Convert from natural log to log2
                'p_value': p_value
            })
        else:
            # Coefficient not found (e.g., due to collinearity)
            results.append({
                'gene': gene,
                'log2FoldChange': np.nan,
                'p_value': np.nan
            })
    results_df = pd.DataFrame(results).set_index('gene')
    return results_df
