#!/usr/bin/env python

# generate a bootstrapped error estimate
# adapted from https://machinelearningmastery.com/calculate-bootstrap-confidence-intervals-machine-learning-results-python/

import numpy as np


def bootstrap_error_estimate(truth, pred, method, random_state=None, alpha=0.95,
                             sample_frac=1.0, iterations=1000, bounds=None):
    """
    Generate a bootstrapped estimate of confidence intervals

    Note the argument order: truth first, then pred, matching the (y_true, y_pred)
    convention used by every scikit-learn metric.

    :param truth: list of experimental values
    :param pred: list of predicted values
    :param method: method to evaluate performance, e.g. matthews_corrcoef
    :param random_state: seed for the resampling, so results are reproducible
    :param alpha: confidence limit (e.g. 0.95 for 95% confidence interval)
    :param sample_frac: fraction of the data to resample on each iteration.  The textbook
        bootstrap resamples n points with replacement, so this defaults to 1.0.  Values
        below 1.0 widen the resulting interval.
    :param bounds: optional (low, high) tuple to clip the interval to the metric's range,
        e.g. (0, 1) for ROC AUC or (-1, 1) for Matthews correlation.  Defaults to no
        clipping, since a metric's range is not knowable from the metric alone.
    :param iterations: number of iterations for resampling
    :return: lower and upper bounds for confidence intervals
    """
    rng = np.random.default_rng(random_state)
    truth = np.asarray(truth)
    pred = np.asarray(pred)
    index_list = np.arange(len(truth))
    num_samples = int(len(index_list) * sample_frac)
    stats = []
    for _ in range(iterations):
        sample_idx = rng.choice(index_list, size=num_samples, replace=True)
        try:
            stats.append(method(truth[sample_idx], pred[sample_idx]))
        except ValueError:
            # A resample can contain only a single class, for which some
            # metrics (e.g. roc_auc_score) are undefined; skip those draws.
            continue
    if len(stats) < 2:
        raise ValueError(
            f"only {len(stats)} of {iterations} bootstrap resamples were usable; "
            "cannot compute a confidence interval"
        )
    p = ((1.0 - alpha) / 2.0) * 100
    lower = np.percentile(stats, p)
    p = (alpha + ((1.0 - alpha) / 2.0)) * 100
    upper = np.percentile(stats, p)
    if bounds is not None:
        lower = max(bounds[0], lower)
        upper = min(bounds[1], upper)
    return lower, upper
