# GPfunctional ------------------------------------------------------------

# manuellement
install.packages("~/M2 DATA ANALYST/STAGE/BCNAM_Data/GPfunctional/",
                 repos = NULL, type="source")
library(GPfunctional)


# FICHIERS ----------------------------------------------------------------

load(file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Global_map.RData")
load(file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/geno_GPfunctional.RData")


# LIBRARY  ----------------------------------------------------------------

library(rrBLUP)
library(regress)


# Temps de calcul (essayé d'optimiser) ------------------------------------

library(profvis)

profvis(expr = mc_GPF_sim_pheno(X = geno, map = map, n_QTL = 100,
                                type_effect = "series", s = 0.9^2))
# 4280ms dont 3520 pour scale center 
# sans scale center on passe à 370 ms


# Tester multi_cross ------------------------------------------------------

# ca marche pour 1 QTL 
d <- mc_GPF_sim_pheno(X = geno, map = map,
                   n_QTL = 10, prop_f = 1, h2f = 1)

var(d$d_y$y_f)

# geno_c1 <- geno[1:100, ]

# Pour 2 mk : ----

# d <- mc_GPF_sim_pheno(X = geno, map = map, 
#                    n_QTL = 2, prop_f = 1, h2f = 1)
# 
# var(d$d_y$y_f)



res <- matrix(NA, nrow = 100, ncol = 4)

for(i in 1:100){
  
  d <- mc_GPF_sim_pheno(X = geno_c2, map = map, 
                        n_QTL = 2, prop_f = 1, h2f = 1)
  res[i, ] <- diag(var(d$d_y))
  
  
}

colnames(res) <- c("Sp", "Sgf", "Sgr", "Se")
res <- res[, -3]
res <- data.frame(res)
boxplot(res)

abline(h = 100, col = "green")
abline(h = 50, col = "blue")

save(res, file= "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/mc_res_1.RData")
# res1 : h2f = 1, propf = 1, 2 QTL, 10 cross

res_mc<-res

# PAS MC 
# res <- matrix(NA, nrow = 100, ncol = 4)
# 
# for(i in 1:100){
#   
#   d <- GPF_sim_pheno(X = geno, map = map, 
#                         n_QTL = 2, prop_f = 1, h2f = 1)
#   res[i, ] <- diag(var(d$d_y))
#   
#   
# }
# 
# colnames(res) <- c("Sp", "Sgf", "Sgr", "Se")
# res <- res[, -3]
# res <- data.frame(res)
# boxplot(res)
# 
# abline(h = 100, col = "green")
# abline(h = 50, col = "blue")
# 
# save(res, file= "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/PASmc_res_1.RData")


## POUR 2 CROSS :
geno_c2 <- geno_c5[1:200, ]
## POUR 5 CROSS :
geno_c5 <- geno[1:500, ]


# check the simulated variance - single cross - 10 mk ----

# geno_c1 <- geno[1:100, ]

# d <- mc_GPF_sim_pheno(X = geno_c1, map = map, n_QTL = 10,
#                    corr_f = 1/sqrt(2))
# 
# var(d$d_y$y_f)
# var(d$d_y$y_r)
# 

res <- matrix(NA, nrow = 50, ncol = 4)
res_GBLUP <- matrix(NA, nrow = 50, ncol = 3)

for(i in 1:20){
  print(i)
  
  d <- mc_GPF_sim_pheno(X = geno_c5, map = map, n_QTL = 100,
                        type_effect = "series", s = 0.9^2)
  
  # d <- mc_GPF_sim_pheno(X = geno_c5, map = map, n_QTL = 100)
  
  res[i, ] <- diag(var(d$d_y))
  
  cr_ind<-substr(x=rownames(geno_c5), start = 1, stop = 4)
  d$d_y$cross <-cr_ind
  
  # variance component estimation
  Kf <- A.mat(X = d$X_f - 1)
  Kr <- A.mat(X = d$X_r - 1)
  
  m <- regress(formula = d$d_y$y_sim ~ as.factor(d$d_y$cross), Vformula = ~ Kf + Kr)
  # m <- regress(formula = d$d_y$y_sim ~ 1, Vformula = ~ Kf + Kr)
  
  res_GBLUP[i, ] <- m$sigma
  
  save(res, file= "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/mc_res_normal_100_series_5cross_fixe.RData")
  save(res_GBLUP, file= "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/mc_res_GBLUP_100_series_5cross_fixe.RData")
  
  
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


# check the simulated variance - single cross - 100 mk ----

res <- matrix(NA, nrow = 20, ncol = 4)
res_GBLUP <- matrix(NA, nrow = 20, ncol = 3)

for(i in 1:20){
  print(i)
  
  d <- mc_GPF_sim_pheno(X = geno, map = map, n_QTL = 100)
  
  res[i, ] <- diag(var(d$d_y))
  
  # variance component estimation
  Kf <- A.mat(X = d$X_f - 1)
  Kr <- A.mat(X = d$X_r - 1)
  
  m <- regress(formula = d$d_y$y_sim ~ 1, Vformula = ~ Kf + Kr)
  res_GBLUP[i, ] <- m$sigma
  
  save(res, file= "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/mc_res_normal_100.RData")
  save(res_GBLUP, file= "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/mc_res_GBLUP_100.RData")
  
  
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


