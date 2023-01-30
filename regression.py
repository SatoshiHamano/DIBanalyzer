import numpy
import scipy.stats


def confband(xd, yd, a, b, x, conf=0.95):
    """
        Calculates the confidence band of the linear regression model at the desired confidence
        level, using analytical methods. The 2sigma confidence interval is 95% sure to contain
        the best-fit regression line. This is not the same as saying it will contain 95% of
        the data points.
        Arguments:
        - conf: desired confidence level, by default 0.95 (2 sigma)
        - xd,yd: data arrays
        - a,b: linear fit parameters as in y=ax+b
        - x: (optional) array with x values to calculate the confidence band. If none is provided, will
        by default generate 100 points in the original x-range of the data.

        Returns:
        Sequence (lcb,ucb,x) with the arrays holding the lower and upper confidence bands
        corresponding to the [input] x array.
        Usage:
        >>> lcb,ucb,x=nemmen.confband(all.kp,all.lg,a,b,conf=0.95)
        calculates the confidence bands for the given input arrays
        >>> pylab.fill_between(x, lcb, ucb, alpha=0.3, facecolor='gray')
        plots a shaded area containing the confidence band
        References:
        1. http://en.wikipedia.org/wiki/Simple_linear_regression, see Section Confidence intervals
        2. http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
        Author: Rodrigo Nemmen
        v1 Dec. 2011
        v2 Jun. 2012: corrected bug in computing dy
        """
    alpha = 1. - conf  # significance
    n = xd.size  # data sample size

    # if x == None: x = numpy.linspace(xd.min(), xd.max(), 100)

    # Predicted values (best-fit model)
    y = a * x + b

    # Auxiliary definitions
    # sd=scatterfit(xd,yd,a,b)	# Scatter of data about the model
    sd = 1. / (n - 2.) * numpy.sum((yd - a * xd - b) ** 2);
    sd = numpy.sqrt(sd)
    sxd = numpy.sum((xd - xd.mean()) ** 2)
    sx = (x - xd.mean()) ** 2  # array

    # Quantile of Student's t distribution for p=1-alpha/2
    q = scipy.stats.t.ppf(1. - alpha / 2., n - 2)

    # Confidence band
    # print q, sd, n, sx, sxd
    dy = q * sd * numpy.sqrt(1. / n + sx / sxd)
    # print x
    ucb = y + dy  # Upper confidence band
    lcb = y - dy  # Lower confidence band

    return lcb, ucb, x


def regressionline(x, y):

    sumx = numpy.sum(x)
    sumy = numpy.sum(y)
    sumxy = numpy.sum(x * y)
    sumxx = numpy.sum(x ** 2)
    sumyy = numpy.sum(y ** 2)
    narray = x.size

    ssxx = sumxx - narray * (sumx / narray) ** 2
    ssyy = sumyy - narray * (sumy / narray) ** 2
    ssxy = sumxy - narray * sumx / narray * sumy / narray
    paa = ssxy / ssxx
    pab = sumy / narray - paa * sumx / narray

    return paa, pab, ccxarray, ccyarray

