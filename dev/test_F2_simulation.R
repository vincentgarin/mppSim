################################
# Test F2 simulation procedure #
################################

library(simcross)
library(mppSim)
library(mppR)
library(LDcorSV)
library(ggplot2)

# sub-functions ----

scale_center_X <- function(X){

  p <- apply(X, 2, FUN = function(x) sum(x)/(2*length(x)))
  P <- rep(1, nrow(X)) %*% (2*t(p))
  Z <- X - P
  den <- sqrt(2 * p * (1 - p))
  Z <- t(t(Z) / den)
  # equivalent to Z %*% diag(1 / den)

  return(Z)

}

# simulation of a 8 crosses NAM population (9 parents) ----

data("geno_par")
rownames(geno_par) <- paste0('P', 1:nrow(geno_par))

data("EUNAM_map")
rownames(map) <- map[, 1]

# NAM crossing scheme with 9 parents
cross_scheme <- cross_scheme_NAM(9)

geno <- sim_mpp_cross(geno_par = geno_par, map = map,
                      cross_scheme = cross_scheme,
                      geno_score_format = "numeric_IBD")

# simulation of the variance of a single QTL ----
X <- geno[1:100, ]
Xs <- scale_center_X(X = X)

var_s <- var_ns <- rep(NA, 100)

# set the genetic variance
Sg <- 50

for(i in 1:100){

  pos_i <- sample(1:ncol(X), size = 1)

  X_f <- X[, pos_i, drop = FALSE]
  Xs_f <- Xs[, pos_i, drop = FALSE]

  # calculate the marker frequency
  p_f <- apply(X_f, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x)))

  # determine the QTL effect
  f <- 1/2
  q_f <- 1 - p_f
  u <- sample(c(1, -1), size = 1)
  Bjf <- u * f *  sqrt((Sg)/(2*p_f*q_f))

  # form the simulated phenotype values
  y_s <- (Xs_f %*% matrix(Bjf, ncol = 1))[, 1]
  y <- (X_f %*% matrix(Bjf, ncol = 1))[, 1]

  var_s[i] <- var(y_s)
  var_ns[i] <- var(y)

}

res <- data.frame(var_ns, var_s)
colMeans(res)
boxplot(res)

abline(h = 12.5, col = "green")
abline(h = 25, col = "blue")

# check the LD pattern ----

# remove the markers that are monomorphic
MAF <- apply(geno, 2, FUN = function(x) sum(x)/(2*length(x)))
rem_mk <- MAF == 0 | MAF == 1
map <- map[!rem_mk, ]
geno <- geno[, !rem_mk]

# select (two) crosses
cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)
geno <- geno[cr_ind %in% unique(cr_ind)[1:2], ]
cr_ind <- cr_ind[cr_ind %in% unique(cr_ind)[1:2]]

# test chr 1
LD_chr1 <- LD_computation(geno = geno, map = map, pop_ind = cr_ind,
                          max_g_dist = 40, chr_sel = 1, prop_mk_comb = 1)

LD_plot(LD_chr1)

LD_decay_01 <- LD_decay(LD_chr1)

# check the potential LD decay when n cross increase ----

geno <- sim_mpp_cross(geno_par = geno_par, map = map,
                      cross_scheme = cross_scheme,
                      geno_score_format = "numeric_IBD")

cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)

LD_dec_thr <- rep(NA, 8)

for(i in 1:8){

  geno_i <- geno[cr_ind %in% unique(cr_ind)[1:i], ]
  cr_ind_i <- cr_ind[cr_ind %in% unique(cr_ind)[1:i]]

  if(length(unique(cr_ind_i)) == 1){cr_ind_i <- NULL }

  # test chr 1
  LD_chr1 <- LD_computation(geno = geno_i, map = map, pop_ind = cr_ind_i,
                            max_g_dist = 40, chr_sel = 1, prop_mk_comb = 1)

  print(LD_plot(LD_chr1))

  LD_dec_thr[i] <- LD_decay(LD_chr1)[1, 2]

}


