# set cell counts, observed externally during experiment
N_treated_t = given, count of BFP-positive at endpoint
N_treated_0 = given, count of BFP-positive at T0
N_untreated_t = given, count of BFP-negative at endpoint
N_untreated_0 = given, count of BFP-negative at T0
t = 30 # given, days
MIN_COUNT_VAL = 1

# Assume dataframe counts contains columns t0, untreated, treated of guide counts
counts = load( count_summary_file )

protect_zero = function( X ){
    # protect against divide by zero; if a guide reduces count to zero from positive 
    # number, don't want to return Infinity
    X[X=0] = MIN_COUNT_VAL
    X
}
# G: growth rate of WT cells
G = 1/t * log2( N_untreated_t / N_untreated_0 )

# K: selective pressure of treatment
K = G - ( 1/t * log2( N_treated_t / N_treated_0 ) )

# gamma, effect of guide in the absence of treatment
count_ratio_untreated = sum(counts$t0) / sum(counts$untreated) 
gamma = 1/(G*t)*( log2( count_ratio_untreated * counts$untreated / protect_zero(counts$t0) ) )

# epsilon, differential enrichment of each guide
count_ratio_treated = sum(counts$untreated) / sum(counts$treated) 
epsilon = log2( count_ratio_treated * counts$treated / protect_zero(counts$untreated) ) 

# rho, resistance a guide confers to selective pressure of treatment
rho = 1 / (K*t) * epsilon
res = cbind( counts, rho, gamma, epsilon )
