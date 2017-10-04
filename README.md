# Arbor 1.0
## by Paul O. Lewis (paul.lewis@uconn.edu)

Program implementing the grape method for estimating the marginal likelihood for phylogenetic models.

The grape method adapts the method described in this paper to phylogenetics:

Wang, Y.-B., M.-H. Chen, L. Kuo, and P. O. Lewis. A new Monte Carlo method for estimating
marginal likelihoods. _Bayesian Analysis_. Published online 28 February 2017.
https://doi.org/10.1214/17-BA1049

## Data File Format

Arbor expects the data file to have the following format:

    iter    tree          log(Like)     log(Prior)   ...
       0       0    -18387.43537599    12.88901141   ...
      10       1    -13867.57430384    11.97246396   ...
      20       2    -13866.47882137    13.11207253   ...
      30       3    -13865.46111627    13.07144537   ...
      40       4    -13867.05433473    9.734989420   ...
      50       5    -13863.91023218    13.82377480   ...
       .       .          .              .            .
       .       .          .              .            .
       .       .          .              .            .

The first four columns are:
1. iter: MCMC iteration
1. tree: index of tree topology
1. log(Like): the natural logarithm of the likelihood
1. log(Prior): the natural logarithm of the joint prior

The remaining columns hold transformed parameter values. Each of these columns must correspond to one of the following types:
* log(edge1), log(edge2), etc.: the natural logarithm of one edge length
* logit(1_pinvar): the natural logarithm of pinvar/(1-pinvar), where pinvar is the proportion of invariable sites
* log(1_gamma_shape): the natural logarithm of the shape parameter of the discrete gamma among-site rate heterogeneity distribution
* log(1_rAG/1_rAC), log(1_rAT/1_rAC), log(1_rCG/1_rAC), log(1_rCT/1_rAC), log(1_rGT/1_rAC): the natural logarithm of the ratio of GTR exchangeabilities to the A-to-C exchangeability
* log(1_freqC/1_freqA), log(1_freqG/1_freqA), log(1_freqT/1_freqA): the natural logarithm of the ratio of each relative nucleotide frequency to the frequency of A

The prefix "1_" in parameter names above refers to the partition subset.
If there were two partition subsets, then columns corresponding to the second
subset would include, e.g., "logit(2_pinvar)", "log(2_gamma_shape)", etc.
Note that Arbor assumes that edge lengths are never unlinked across partition subsets.
Also, the "1_" prefix must be used even for unpartitioned data.
