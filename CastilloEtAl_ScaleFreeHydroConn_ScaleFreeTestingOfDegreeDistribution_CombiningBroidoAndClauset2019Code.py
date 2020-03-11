'''
This script is provided as part of the methods section for the manuscript entitled "Scale-free
structure of surface-water connectivity within a lowland river floodplain" by Cesar R Castillo,
İnci Güneralp, Billy Hales, and Burak Güneralp. The corresponding author is Cesar R Castillo
and can be contacted at castillocesar@tamu.edu.  

This script uses the degree sequences created by the 
CastilloEtAl_ScaleFreeHydroConn_CreateDegreeSequenceFromDegreeDistribution_ForBroidoAndClauset2019Code
code and subjects it to code provided by Broido and Clauset in their 2019 paper published in
Nature Communications under the title "Scale-free networks are rare". This script combines all their
code and reads in the degree sequence and conducts an empirical analysis on the degree distribution
that is used to determine how scale-free a network is.

The original Broido and Clauset 2019 code can be found at the URL below
https://github.com/adbroido/SFAnalysis/tree/master/code

This code was designed to be used in Python 3.6 with the most up-to-date version of the libraries
listed below
Python libraries: pandas, numpy, scipy, time, mpmath, norm, chi2, and os
'''



'''
code for importing the a degree sequence of a graph/network
the original code itself can be downloaded from the URL below
https://github.com/adbroido/SFAnalysis/blob/master/code/importfiles.py
'''

import pandas as pd
import numpy as np


def readdata(fp):
    """
    Reads in a datafile.
    Input:
        fp                      string, filepath to csv file (degree sequence).
    Output:
        data                    ndarray, ndim = 1, dtype = integer. Repeats each
                                xval as many times as indicated by counts
    """

    df = pd.read_csv(fp)
    data = np.concatenate([np.array([df.xvalue[i] for repeat in range(df.counts[i])]) for i in range(len(df))])
    return data
# end of function

##########################################################################################################################

"""
code for fitting heavy-taile distributions to node degrees from a 
supposed scale free graph 
code is from the URL below
https://github.com/adbroido/SFAnalysis/blob/master/code/fit.py
"""

import numpy as np
import scipy.optimize as op
#import integration_constants as ic
import scipy.special as sp
import time

"""
Contains functions used in fitting the power-law, exponential, log-normal,
Weibull (stretched exponential), and power-law with exponential cutoff, as well
as plpval() to find a p-value for the power-law fit. All should be called
directly. All distributions are discrete.
"""


def pl(x):
    """ Fits a tail-conditional power-law to a data set. This implements brute
    force optimization (grid search) instead of using a built in optimizer. The
    grid on alpha runs from alstart to alstart+shift. This is based on Aaron's
    plfit.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).
    Input:
        x           ndarray, ndim = 1, dtype = integer
    Output:
        alpha        float, exponent on x, must be > 1
        xmin         int, starting point for power law tail, must be >= 1
        ntail        int, number of datapoints above (and including) xmin
        L            float, log likelihood of the returned fit
        ks           float, goodness of fit statistic (Kolmogorov-Smirnov)
    """
    # find the and sort unique possible xmin values
    xminV = np.trim_zeros(np.unique(x))
    # initialize array of the fits for every xmin
    fitV = np.zeros([len(xminV),2])

    start_time = time.time()
    # initialize vector of constants
    # where the xmins start
    xminprev = min(xminV) - 1
    # initialize array of possible alpha values
    alstart = 1.01
    shift = 9.50
    alphaV = np.arange(alstart,alstart+shift,0.01)
    zetaV = sp.zeta(alphaV)
    constV = zetaV
    # shift up to start at the smallest xmin
    for j in range(xminprev):
        constV += -(1+j)**(-alphaV)

    # loop over the xmin values at find the best fit at each
    for i in range(len(xminV)):
        xmin = xminV[i]
        xtail = x[x>=xmin]
        ntail = len(xtail)
        # optimize over alpha
        # find the corresponding array of conditional log likelihoods
        Ls = -alphaV*np.sum(np.log(xtail)) - ntail*np.log(constV)
        # pull out the location of the best fit alpha
        aind = Ls.argmax()
        # find what alpha value is at this index
        alpha = alphaV[aind]
        # compute the KS statistic
        # theoretical cdf
        cdf = np.cumsum(range(np.min(xtail), np.max(xtail)+1)**(-alpha)/constV[aind])
        #  binned data
        xhist = np.histogram(xtail,range(np.min(xtail), np.max(xtail)+2))
        # empirical cdf
        edf = np.cumsum(xhist[0])/float(ntail)
        # KS stat
        ks = np.max(np.abs(cdf-edf))
        # add this KS stat and alpha to the array of fits
        fitV[i] = np.array([ks, alpha])
        # update the constants
        for j in range(xmin-xminprev):
            constV += -(xmin+j)**(-alphaV)
        xminprev = xmin


    # pull out the index of the smallest KS stat
    ksind = fitV[:,0].argmin()
    ks = fitV[ksind,0]
    # find the corresponding xmin
    xmin = xminV[ksind]
    # find the corresponding alpha
    alpha = fitV[ksind,1]
    # evaluate the likelihood here
    xtail = x[x>=xmin]
    ntail = len(xtail)
    start_time = time.time()
    const = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    L = -alpha * np.sum(np.log(xtail)) - ntail*np.log(const)
    # print "-------%s seconds -----------" %(time.time()-start_time)
    # print "alpha = %s" %alpha
    # print "xmin = %s" %xmin
    return [alpha,xmin, ntail, L, ks]
# end of function

def plpval(x, alpha, xmin, gof):
    """
    Finds p-value for the power-law fit using a KS test. This is based on
    Aaron's plpva.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).
    Input:
        x            ndarray, ndim = 1, dtype = integer
        alpha        float, exponent on x, must be > 1
        xmin         int, starting point for power law tail, must be >= 1
        gof           float, goodness of fit statistic (Kolmogorov-Smirnov)
    Output:
        p            p-value of the returned fit (reject PL hypothesis for p<0.1)
    """
    # set desired precision level in p-value
    eps = 0.01
    #num_resamps = int(np.ceil((1./4)*eps**(-2)))
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    xmax = np.max(x)
    tailinds = x>=xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail)/n
    mmax = 20*xmax
    # set the tail of the pdf
    #const_tail = ic.plconst(np.array(alpha),xmin)
    const_tail = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    pdf_tail = np.arange(xmin,mmax+1)**(-alpha)/const_tail # i.e.; end at mmax
    # pad this with zeros (rewrite if we don't need to do this)
    pdf = np.zeros(mmax+1)
    pdf[xmin:] = pdf_tail
    # clean up in case this is a huge array
    del pdf_tail
    # set the cdf. rows are x-val or cdf(xval). So cdf(x=10) is cdf[1,10]
    cdf = np.array( [ np.arange(mmax+1), np.cumsum(pdf) ] )
    # tack on a last entry
    cdf = np.concatenate( (cdf , np.array([[mmax+1,1]]).T) , axis = 1 )

    # semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # non-parametric bootstrap from the head of x
        # count how many of n random numbers are in the head, based on the probability of being in the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n)>ptail)
        headinds = np.array([np.floor(nhead*np.random.rand(nnewhead))],dtype=int)
        newhead = xhead[headinds][0]
        nnewtail  = n-nnewhead

        # parametric bootstrap for the powerlaw tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype = int)
        indrtail = 0
        indnewtail = 0
        for xval in range(xmin, mmax+2):
            while (indrtail < len(rtail)) and (rtail[indrtail] <= cdf[1, xval]):
                indrtail += 1
            newtail[indnewtail:indrtail] = xval
            indnewtail = indrtail
            if indnewtail > nnewtail:
                break
        # combine into new sample
        newx = np.concatenate((newhead, newtail))
        if (newx == np.zeros_like(newx)).all():
            import pdb; pdb.set_trace()
        # fit this new sample
        [newalpha, newxmin, newntail, newLpl, newgof] = pl(newx)
        # print where we are
        current_p = np.sum(bootstraps[0:resamp_ind]>=gof)/(float(resamp_ind+1))
        # print "[%s]    p = %f" %(resamp_ind, current_p)
        # store gof stat
        bootstraps[resamp_ind] = newgof
        # if it's taking forever and we can end, do it
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps/20.:
                if current_p<0.05 or current_p>0.5:
                    print("current p = %s   elapsed time = %s" %(current_p, time.time()-starttime))
                    return current_p
    p = np.sum(bootstraps>=gof)/float(num_resamps)
    print("p = %s   elapsed time = %s" %(p, time.time()-starttime))
    return p
# end of function

def exp(x):
    """
    Fits a tail-conditional exponential to a data set. The data is assumed
    to begin at xmin. The logpdf is what is calculated and returned, as this is
    more relevant for likelihood calculations.
    Input:
        x            ndarray, ndim = 1, dtype = integer
    Output:
        lam          float, exponential rate, must be > 0
        LV           ndarray, pointwise log likelihood of the returned fit
        convstatus   Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def logpdf(x,lam):
        result = np.log(1-np.exp(-lam))+lam*xmin - lam*x
        return result
    if len(np.unique(x))<2:
        # this means every value is equal to xmin
        # return dummy answers and say we don't converge
        lam = 0
        LV =0
        convstatus = False
    else:
        # Moment based estimate for optimzation
        lam0 = np.log(1+float(ntail)/np.sum(x-xmin))
        # define negative log likelihood, the function we wish to minimize
        negloglike = lambda lam: -np.sum(logpdf(x,lam))
        tol = 1E-9
        res = op.minimize(negloglike,lam0, bounds=[(tol,None)],method='L-BFGS-B')
        lam = np.asscalar(res.x)
        convstatus = res.success
        LV = logpdf(x,lam)
    return [lam, LV, convstatus]
# end of function

def ln(x):
    """
    Fits a tail-conditional log normal distribution to a data set.
    The data is assumed to begin at xmin. The logpdf is what is calculated and
    returned, as this is more relevant for likelihood calculations.
    Discretization is done by binning the continuous distrbution
    (see text for details)
    Input:
        x               ndarray, ndim = 1, dtype = integer
    Output:
        theta           ndarray, [mu, sigma] where mu is a float, the mean of the
                            distribution, unbounded (though we impose a bound), and sigma is a
                            float, the standard deviation, and must be > 0
        LV              ndarray, pointwise log likelihood of the returned fit
        convstatus      Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def logpdf(x, mu, sigma):
        xmin = np.min(x)
        F = lambda x: (sp.erfc((np.log(x)-mu)/(np.sqrt(2)*sigma)))/2
        g = lambda x: F(x)- F(x+1)
        h = -np.log(F(xmin))+np.log(g(x))
        return h
    # initial estimates
    mu0 = 0
    sigma0 = 1
    theta0 = np.array([mu0, sigma0])
    n = len(x)
    # optimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-1
    bnds=[(-n/5,None),(tol,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds, method='L-BFGS-B')
    theta = res.x
    convstatus = res.success
    LV = logpdf(x,theta[0], theta[1])
    return [theta, LV, convstatus]
# end of function

'''
def plwc(x, alpha0=None):
    """
    Fits a tail-conditional power-law with exponential cutoff to a data set.
    The data is assumed to begin at xmin. The logpdf is what is calculated and
    returned, as this is more relevant for likelihood calculations.
    Input:
        x           ndarray, ndim = 1, dtype = integer
        alpha0      float, power-law exponent (optional input)
    Output:
        alpha        float, exponent on x, must be > -1
        lam          float, exponential rate, must be > 0
        LV           ndarray, pointwise log likelihood of the returned fit
        convstatus   Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def logpdf(x,alpha, lam):
        xmin = np.min(x)
        C = ic.plwcconst(alpha,lam, xmin)
        result = -np.log(C) - alpha*np.log(x) - lam*x
        return result
    # Estimates for optimzation
    if alpha0 is None:
        alpha0 = pl(x)[0]
    lam0 = exp(x)[0]
    theta0 = np.array([alpha0,lam0])
    # define negative log likelihood, the function we wish to minimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0], theta[1]))
    tol = 1E-5
    bnds=[(-1+tol,None),(tol,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds)
    # res = op.minimize(negloglike,theta0, method='Nelder-Mead')
    theta = res.x
    convstatus = res.success
    alpha = theta[0]
    lam = theta[1]
    LV = logpdf(x,alpha, lam)
    return [alpha, lam, LV, convstatus]
# end of function
'''


def strexp(x):
    """
    Fits a tail-conditional stretched exponential distribution to a data set.
    The data is assumed to begin at xmin. The logpdf is what is calculated and
    returned, as this is more relevant for likelihood calculations.
    Discretization is done by binning the continuous distrbution
    (see text for details)
    Input:
        x           ndarray, ndim = 1, dtype = integer
    Output:
        theta           ndarray, [a,b], dtype=float. a (0<a<1) is the mean of the
                            distribution and b is the standard deviation (b > 0)
        LV              ndarray, pointwise log likelihood of the returned fit
        convstatus      Boolean, True if the fit converged, false if not
    """
    xmin = np.min(x)
    ntail = len(x)
    # define log pdf
    def initialguessweib(x):
        """
        The follow is not an efficient estimator of the shape, but serves to start
        the approximation process off.  (See Johnson and Katz, ch. 20, section 6.3,
        "Estimators Based on Distribution of log X".) This is taken from the R
        package by Cosma Shalizi at http://tuvalu.santafe.edu/~aaronc/powerlaws/).
        """
        xmin = np.min(x)
        n = len(x)
        shape = (np.sqrt(6)/np.pi)*np.std(np.log(x))
        scale = (np.sum(x**shape)/n)**(1/shape)
        return np.array([shape,scale])

    # define log pdf, for use in likelihood calculation
    def logpdf(x, a, b):
        xmin = np.min(x)
        F = lambda x: np.exp(-(x/b)**a)
        g = lambda x: F(x)-F(x+1)
        h = -np.log(F(xmin))+np.log(g(x))
        return h
    # initial estimates
    # initial estimates
    theta0 = initialguessweib(x)
    # optimize
    negloglike = lambda theta: -np.sum(logpdf(x,theta[0],theta[1]))
    tol = 1E-5
    bnds=[(tol,1),(0.01,None)]
    res = op.minimize(negloglike, theta0, bounds=bnds, method='L-BFGS-B')
    theta = res.x
    convstatus = res.success
    LV = logpdf(x,theta[0], theta[1])
    return [theta, LV, convstatus]
# end of function

#####################################################################################################################

'''
code for calculating the normalization constants in heavy-tailed distributions
the code itself is from the URL below
https://github.com/adbroido/SFAnalysis/blob/master/code/integration_constants.py
'''

import mpmath as mp
import numpy as np

""" Contains functions to compute the normalization constants for both power-law
and power-law with expoenential cutoff. Both functions are meant to be called
directly.
Note that plwcconstant returns a constant to divide by, while plconst returns
a constant to multiply.
"""


def plwcconst(alpha, lam, xmin):
    """ Computes the normalization constant on the discrete power law;
    i.e., computes C so that
        1/C * sum from xmin to infinity of ( x^(-alpha) * e^(-lam x)] ) = 1
    The formula below is obtained by noting that the sum above is equal to
        sum from xmin to infinity of ( x^(-alpha) * e^(-lam x)] ) =
        e^(-xmin * lam) * HurwitzLerchPhi(e^(-lam), alpha, xmin)
    where HurwitzLerchPhi is the Lerch Phi function:
        Phi(z,s,a) = sum from 0 to infinity of ( z^k / (a+k)^s )
    (as defined in sympy docs). If one is disinclined to note this fact
    by deriving it, one is encouraged to look at it in Mathematica and
    blindly trust the result.
    Note: The standard notation for the Lerch Phi function above uses "a" as
          a parameter. This does not correspond to the alpha we pass to the
          function.
    Inputs:
        alpha                  float, exponent on x, must be > -1
        lam                    float, exponential cutoff, must be > 0
        xmin                   int, starting point for sum, must be >= 1
    Outputs:
        C                      float, normalization constant
    """
    mp.mp.dps = 40
    result = mp.exp(-xmin * lam) * mp.lerchphi(mp.exp(-lam), alpha, xmin)
    C = float(result)
    return C
# end of function

def plconst(alpha, xmin):
    """ Computes the normalization constant on the discrete power law;
    i.e., computes C so that
        C * sum from xmin to infinity of x^(-alpha) = 1
    The formula below is obtained by noting that the sum above is equal to
        sum from xmin to infinity of x^(-aplha) = HurwitzZeta(alpha, xmin)
    where HurwitzZeta is defined by
        zeta(s,a) = sum from 0 to infinity of 1 / (a+k)^s
    (as defined in sympy docs). If one is disinclined to note this fact
    by deriving it, one is encouraged to look at it in Mathematica and
    blindly trust the result.
    Note: The standard notation for the zeta function above uses "a" as
          a parameter. This does not correspond to the alpha we pass to the
          function. (Our alpha is s in the notation above).
    Inputs:
        alpha                  array, shape=(1,) exponent on x, must be > 1
                               must be passed as array for op.minimize()
        xmin                   int, starting point for sum, must be >= 1
    Outputs:
        C                      float, normalization constant
    """
    total = mp.zeta(np.asscalar(alpha),1) # op.minimize passes array
    lowertail = np.sum(np.asarray(range(1,xmin)**(-alpha)))
    result = total-lowertail
    C = 1./(result)
    return float(C)
# end of function


# example usage:
if __name__ == '__main__':
    alpha = np.array([1.2])
    lam = 2.0
    xmin = 3

    const1 = plwcconst(alpha, lam, xmin)
    const2 = plconst(alpha, xmin)
# end of code block

##########################################################################################################################

'''
this code is used to test a power law fit with other heavy-tailed distributions
this code was copied from the URL below
https://github.com/adbroido/SFAnalysis/blob/master/code/lrt.py
'''

import numpy as np
import scipy.optimize as op
import scipy.special as sp
#import integration_constants as ic
import fit
from scipy.stats import norm, chi2


""" Contains functions used in likelihood ratio tests, comparing the power-law
fit with alternative distributions. All can be called directly, but nested() and
nonnested() are the only we use.
"""

'''
def pllogpdf(x,alpha):
    """ Point-wise log-pdf of the power-law distribution. For use in computing
    point-wise log-likelihood ratios.
    Input:
        x               ndarray, data to be fit
        alpha           float, power-law exponent, comes from best-fit
    Output:
        logpdf          float, point-wise log pdf
    """
    logpdf = np.log(ic.plconst(np.array([alpha]), np.min(x))) -alpha*np.log(x)
    return logpdf
# end of function
'''

def vuong(LplV, LaltV):
    """ Vuong test for log-likelihood ratio tests. Computes the ratio and an
    associated p-value.
    Input:
        LplV            ndarray, power-law pointwise log-likelihood
        LaltV           ndarray, alternative pointwise log-likelihood
    Output:
        R               float, likelihood ratio
        p2              float, 2-tail p-value
        normR           float, normalized ratio (so we can assume standard
                            normal for calculating p-val)
    """
    logratioV = LplV - LaltV
    n = len(logratioV)
    R = np.sum(logratioV)
    # standard deviation of normal dist
    sigma = np.std(logratioV)
    normR = (1/np.sqrt(n))*R/sigma
    # one-sided p-value
    p1 = norm.cdf(normR)
    if p1 > 0.5:
        p1 = 1-p1
    # 2-sided p-value
    p2 = 2*p1
    return R, p2, normR
# end of function

def decide(logratio, p, decisionthresh):
    """ Takes p-value and log-likelihood ratio and makes a decision whether the
    test is conclusive. If so, it indicates whether the power-law or alternative
    is favored.
    Input:
        logratio            float, normalized log-likelihood ratio
        p                   float, associated p-value
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1
    Output:
        d                   int, decision from the LRT
                                Decisions are:
                                    1   -  power law better
                                    0   -  inconclusive
                                   -1   -  alternative dist better
    """
    if p <= decisionthresh:
        if logratio == 0:
            d = 0
        elif logratio > 0:
            # PL better
            d = 1
        else:
            # alt better
            d = -1
    else:
        # inconclusive
        d = 0
    return d
# end of function

def decidenested(logratio, p, decisionthresh):
    """ Takes p-value and log-likelihood ratio and makes a decision whether the
    test is conclusive. If so, it indicates whether the power-law or alternative
    is favored. This is for nested decision, which cannot come out in favor of
    the power-law. The only nested alternative distribution we use is power-law
    with exponential cutoff.
    Input:
        logratio            float, normalized log-likelihood ratio
        p                   float, associated p-value
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1
    Output:
        d                   int, decision from the LRT
                                Decisions are:
                                    1   -  power law better
                                    0   -  inconclusive
                                   -1   -  alternative dist better
    """
    if p <= decisionthresh and logratio < 0:
            # alt better
            d = -1
    else:
        # inconclusive
        d = 0
    return d
# end of function

def exp(x, LplV, decisionthresh):
    """
    Perform likelihood ratio test for exponetial distribution. First fits an
    exponential distribution to the data. A Vuong statistic is calculated from
    the likelihoods of the exponential fit and the power law fit.
    Input:
        x                   ndarray, data set to be fit (has been truncated to start
                                at xmin)
        LplV                ndarray, pointwise likelihood values for power-law fit
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1
    Output:
        dexp                int, decision about exponential distribution
    """
    # perform lrt: Log-likelihood ratio between discrete power law and
    # exponential distribution. This is done pointwise so that we can use
    # Vuong's statistic to estimate the variance in the ratio
    [lam, LexpV, convstatus] = fit.exp(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LexpV)
        # check if statistically significant
        dexp = decide(normR, p, decisionthresh)
    else:
        dexp = 2
    return dexp
# end of function

def ln(x,LplV, decisionthresh):
    """
    Perform likelihood ratio test for log normal distribution. First
    fits a log normal distribution to the data. A Vuong statistic is
    calculated from the likelihoods of the exponential fit and the power law
    fit.
    Input:
        x                   ndarray, data set to be fit
        LplV                ndarray, pointwise loglikelihood for power-law fit
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1
    Output:
        dln                 int, decision about log-normal distribution
    """
    [theta,LlnV, convstatus] = fit.ln(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LlnV)
        # check if statistically significant
        dln = decide(normR, p, decisionthresh)
    else:
        dln = 2
    return dln
# end of function

def strexp(x,LplV, decisionthresh):
    """
    Perform likelihood ratio test for stretched exponetial (Weibull)
    distribution. First fits a stretched exponential distribution to the data,
    then calculates a Vuong statistic from the likelihoods of the exponential
    fit and the power law fit.
    Input:
        x                   ndarray, data set to be fit
        LplV                ndarray, pointwise loglikelihood for power-law fit
        decisionthresh      float, threshold for rejecting.
                                Default in paper is decisionthresh = 0.1
    Output:
        dstrexp             int, decision about exponential distribution
    """
    [theta, LstrexpV, convstatus] = fit.strexp(x)
    if convstatus == True:
        R, p, normR = vuong(LplV, LstrexpV)
        # check if statistically significant
        dstrexp = decide(normR, p, decisionthresh)
    else:
        dstrexp = 2
    return dstrexp
# end of function

def nested(x, alpha, decisionthresh=0.1):
    """
    Perform likelihood ratio tests for alternative distributions that are in
    the power law family.
    Input:
        x                   ndarray, data set to be fit. Assumed to only have values above
                                the xmin value from the power law fit
        alpha               float, best fit power-law parameter
        decisionthresh      float, threshold for rejecting null hypothesis
                                Default is 0.1
    Output:
        dplwc    int, decision about power law with exponential cutoff
                    Decisions are   1  -   power law better
                                    0  -   inconclusive
                                   -1  -   alternative dist better
    """
    LplV = pllogpdf(x,alpha)
    Lpl = np.sum(LplV)
    # compare plwc
    [alpha, lam, LplwcV, convstatus] = fit.plwc(x, alpha)
    if convstatus == True:
        Lplwc = np.sum(LplwcV)
        R = Lpl-Lplwc
        p = 1-chi2.cdf(-2*R, df=1)
        dplwc = decidenested(R, p, 0.1)
    else:
        dplwc = 2
    return dplwc
# end of function

def nonnested(x, alpha, decisionthresh=0.1):
    """
    Perform likelihood ratio tests for alternative distributions that are not
    in the power law family.
    Input:
        x                   ndarray, data set to be fit. Assumed to only have
                                values above the best fit xmin value (from PL)
        alpha               float, best fit power-law parameter
        decisionthresh      float, threshold for rejecting null hypothesis
                                        Default is 0.1
    Output:
        dexp                int, decision about exponential distribution
        dln                 int, decision about log normal distribution
        dstrexp             int, decision about stretched exp distribution
                                Decisions are   1  -   power law better
                                                0  -   inconclusive
                                               -1  -   alternative dist better
    """
    LplV = pllogpdf(x,alpha)
    # compare exponential
    dexp = exp(x,LplV, decisionthresh)
    # compare log normal
    dln = ln(x,LplV, decisionthresh)
    # compare stretched exponential
    dstrexp = strexp(x,LplV, decisionthresh)

    return [dexp, dln, dstrexp]
# end of function

#################################################################################################################

"""
code is for testing if the graphs from the inundation analysis are scale free
the statistical code is from the the paper listed below
Broido, A. D., and A. Clauset (2019), Scale-free networks are rare, Nature Communications, 10(1).
"""

import os

# defining the input and output filepaths
# general directory where the inputs and outputs will be
# held
gen_path = 'DIRECTORY'
# subdirectory that contains the set of input degree sequences
in_degseq_path = os.path.join(patch_path,
                              'SUBDIRECTORY')

#########################################################################################################################

"""
creating a list of numpy arrays that hold the sets of degree
sequences that will be processed using the remaining code
"""
import os

files = os.listdir(in_degseq_path)
degseq_list = []
for f in files:
    dat = readdata(os.path.join(in_degseq_path,
                                f))
    degseq_list.append(dat.astype(int))

# end of loop

###########################################################################################################################

"""
creating a list of values that pertain to the degree sequences
this particular code block fits a power law to the data
"""

# creating a vector that depicts the stage height at gaging station
# here we use values that span the historical range for the 08189500
# USGS gaging station on Mission River 
# change as needed for one's own respective analysis
stage = np.linspace(start=float(files[0].split('_')[2][:4])/100,
                    stop=float(files[-2].split('_')[2][:4])/100,
                    num=len(files[:-1]))
stage = np.append(stage,
                  11.63)
                 
plfit_list = []
for i in range(0, len(degseq_list)):
    plfit = pl(degseq_list[i])
    plfit_p = plpval(x=degseq_list[i],
                     alpha=plfit[0],
                     xmin=plfit[1],
                     gof=plfit[4])
    plfit = np.append(plfit, [plfit_p, stage[i]])
    #plfit = [stage[i]] + [plfit]
    plfit_list.append(plfit)
# end of for loop


"""
creating a list of values that pertain to the exponential fit
of the same data used in the power law fitting
"""

expfit_list = []
for i in range(0, len(degseq_list)):
    ar = degseq_list[i]
    ar = ar[ar >= plfit_list[i][1]]
    expfit = exp(ar)
    expfit = [expfit] + [stage[i]]
    expfit_list.append(expfit)
# end of for loop


"""
creating a list of values that pertain to the lognormal fit
of the same data used in the power law fitting
"""

lnfit_list = []
for i in range(0, len(degseq_list)):
    ar = degseq_list[i]
    ar = ar[ar >= plfit_list[i][1]]
    lnfit = exp(ar)
    lnfit = [lnfit] + [stage[i]]
    lnfit_list.append(lnfit)
# end of for loop


"""
creating a list of values that pertain to the tail-conditional
stretched exponential fit of the same data used in the
power law fitting
"""

strexpfit_list = []
for i in range(0, len(degseq_list)):
    ar = degseq_list[i]
    ar = ar[ar >= plfit_list[i][1]]
    strexpfit = strexp(ar)
    strexpfit = [strexpfit] + [stage[i]]
    strexpfit_list.append(strexpfit)
# end of for loop



