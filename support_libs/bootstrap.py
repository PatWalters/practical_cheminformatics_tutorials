#!/usr/bin/env python

# generate a bootstrapped error estimate
# adapted from https://machinelearningmastery.com/calculate-bootstrap-confidence-intervals-machine-learning-results-python/

from __future__ import print_function
import numpy as np
from sklearn.utils import resample


def bootstrap_error_estimate(pred, truth, method, method_name="", alpha=0.95, sample_frac=0.5, iterations=1000):
    """
    Generate a bootstrapped estimate of confidence intervals
    :param pred: list of predicted values
    :param truth: list of experimental values
    :param method: method to evaluate performance, e.g. matthews_corrcoef
    :param method_name: name of the method for the progress bar
    :param alpha: confidence limit (e.g. 0.95 for 95% confidence interval)
    :param sample_frac: fraction to resample for bootstrap confidence interval
    :param iterations: number of iterations for resampling
    :return: lower and upper bounds for confidence intervals
    """
    index_list = range(0, len(pred))
    num_samples = int(len(index_list) * sample_frac)
    stats = []
    for _ in range(0, iterations):
        sample_idx = resample(index_list, n_samples=num_samples)
        pred_sample = [pred[x] for x in sample_idx]
        truth_sample = [truth[x] for x in sample_idx]
        stats.append(method(pred_sample, truth_sample))
    p = ((1.0 - alpha) / 2.0) * 100
    lower = max(0.0, np.percentile(stats, p))
    p = (alpha + ((1.0 - alpha) / 2.0)) * 100
    upper = min(1.0, np.percentile(stats, p))
    return lower, upper

