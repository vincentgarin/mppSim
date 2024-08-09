library(GPfunctional)

data("geno_par")
rownames(geno_par) <- paste0('P', 1:nrow(geno_par))

data("EUNAM_map")
rownames(map) <- map[, 1]

# NAM crossing scheme with 9 parents
cross_scheme <- cross_scheme_NAM(9)

geno <- sim_mpp_cross(geno_par = geno_par, map = map,
                      cross_scheme = cross_scheme,
                      n_ind_cr = rep(100, nrow(cross_scheme)))

# check the simulated variance - single cross - 1 mk ----

geno_c1 <- geno$geno_num_IBD[1:100, ]
colnames(map)[1] <- "mk.names"

d <- GPF_sim_pheno(X = geno_c1, map = map, n_QTL = 1, prop_f = 1, h2f = 1,
                   corr_f = 1/sqrt(2))

var(d$d_y$y_f)

res <- matrix(NA, nrow = 100, ncol = 4)

for(i in 1:100){

  d <- GPF_sim_pheno(X = geno_c1, map = map, n_QTL = 1, prop_f = 1, h2f = 1,
                     corr_f = 1/sqrt(2))
  res[i, ] <- diag(var(d$d_y))


}

colnames(res) <- c("Sp", "Sgf", "Sgr", "Se")
res <- res[, -3]
res <- data.frame(res)
boxplot(res)

abline(h = 100, col = "green")
abline(h = 50, col = "blue")

# check the simulated variance - single cross - 10 mk ----

geno_c1 <- geno$geno_num_IBD[1:100, ]

d <- GPF_sim_pheno(X = geno_c1, map = map, n_QTL = 10,
                   corr_f = 1/sqrt(2))

var(d$d_y$y_f)
var(d$d_y$y_r)

res <- matrix(NA, nrow = 100, ncol = 4)
res_GBLUP <- matrix(NA, nrow = 100, ncol = 3)

for(i in 1:100){

  d <- GPF_sim_pheno(X = geno_c1, map = map, n_QTL = 10,
                     corr_f = 1/sqrt(2))
  res[i, ] <- diag(var(d$d_y))

  # variance component estimation
  Kf <- A.mat(X = d$X_f - 1)
  Kr <- A.mat(X = d$X_r - 1)

  m <- regress(formula = d$d_y$y_sim ~ 1, Vformula = ~ Kf + Kr)
  res_GBLUP[i, ] <- m$sigma

}

# data variance components
colnames(res) <- c("Sp", "Sgf", "Sgr", "Se")
res <- data.frame(res)
boxplot(res)

abline(h = 100, col = "green")
abline(h = 50, col = "blue")
abline(h = 25, col = "red")

# BLUP estimation
colnames(res_GBLUP) <- c("Sgf", "Sgr", "Se")
res_GBLUP <- data.frame(res_GBLUP)

boxplot(res_GBLUP)
abline(h = 25, col = "red")
abline(h = 50, col = "blue")

# variance component estimation with IBS data ----

res <- matrix(NA, nrow = 100, ncol = 4)
res_GBLUP <- matrix(NA, nrow = 100, ncol = 3)

X_IBS <- geno$geno_num_IBS[1:100, ]

for(i in 1:100){

  d <- GPF_sim_pheno(X = geno_c1, map = map, n_QTL = 10,
                     corr_f = 1/sqrt(2))
  res[i, ] <- diag(var(d$d_y))

  # variance component estimation
  X_IBS_f <- X_IBS[, d$mk_sel_f[, 1]]
  X_IBS_r <- X_IBS[, d$mk_sel_r[, 1]]

  Kf <- A.mat(X = X_IBS_f - 1)
  Kr <- A.mat(X = X_IBS_r - 1)

  m <- tryCatch(regress(formula = d$d_y$y_sim ~ 1, Vformula = ~ Kf + Kr),
                error = function(e) NULL)
  if(!is.null(m)){
    res_GBLUP[i, ] <- m$sigma
  }

}

# data variance components
colnames(res) <- c("Sp", "Sgf", "Sgr", "Se")
res <- data.frame(res)
boxplot(res)

abline(h = 100, col = "green")
abline(h = 50, col = "blue")
abline(h = 25, col = "red")

# BLUP estimation
colnames(res_GBLUP) <- c("Sgf", "Sgr", "Se")
res_GBLUP <- data.frame(res_GBLUP)

boxplot(res_GBLUP)
abline(h = 25, col = "red")
abline(h = 50, col = "blue")
