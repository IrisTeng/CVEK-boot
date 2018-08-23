# row 1-63 are linear case, 66-135 quadratic, 138-207 cubic and 211-280 rbf
# cubic and rbf just slightly changed from quadratic (kern and data-generating funs)
#
# take quadratic as example: 
# kernel matches data-generating mechanism (row 80-87);
# after row 86 and 91, ncol of X becomes 3, while it should be 5 (including intercept)
# when I try to fit kernel regression, first I generated kernel like row 107-108,
# but the result was different and test error was pretty large,
# so I changed them into row 109-110 (Is it correct to do so?)
#
# If we don't standardize X and K in CVEK, the estimated projection matrix A is the
# same as the one in row 126. However, if we use svd to compute K, we might have new
# trouble. 
