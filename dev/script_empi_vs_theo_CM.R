# Pour 2 croisements  ------------------------------------------------------

## Exemple pour deux croisements ----

N <- 309 # Nombre d'individus

# Fréquences et effets pour chaque croisement
p1 <- 0.884375; q1 <- 1 - p1; b <- 1
p2 <- 0.614094; q2 <- 1 - p2; c <- 2

# Probabilités d'être dans chaque croisement
phi1 <- 0.5177994
phi2 <- 1 - phi1

# Simuler les génotypes et phénotypes pour chaque croisement
genotypes1 <- sample(c("AA", "AB", "BB"), N * phi1, 
                     replace = TRUE, prob = c(p1^2, 2*p1*q1, q1^2))
genotypes2 <- sample(c("AA", "AC", "CC"), N * phi2, 
                     replace = TRUE, prob = c(p2^2, 2*p2*q2, q2^2))
# proba sous HW

phenotypes1 <- ifelse(genotypes1 == "AA", 0, 
                      ifelse(genotypes1 == "BB", b, b/2))
phenotypes2 <- ifelse(genotypes2 == "AA", 0, 
                      ifelse(genotypes2 == "CC", c, c/2))

# Calculer les variances empiriques pour chaque croisement
# var_empirique1 <- var(phenotypes1)
# var_empirique2 <- var(phenotypes2)

# Calculer la variance globale empirique
phenotypes_total <- c(phenotypes1, phenotypes2)
variance_empirique_totale <- var(phenotypes_total)

# Calculer les variances théoriques pour chaque croisement
# var_theorique1 <- 0.5 * b^2 * q1 * p1
# var_theorique2 <- 0.5 * c^2 * q2 * p2

# Calculer la variance théorique globale
C1 <- (0.5 * p1 * phi1 * q1) + (q1^2 * (phi1 * (1 - phi1)))
C2 <- (0.5 * p2 * phi2 * q2) + (q2^2 * (phi2 * (1 - phi2)))
C3 <- q1 * q2 * phi1 * phi2

variance_theorique_totale <- (b^2 * C1) + (c^2 * C2) - (2 * b * c * C3)

# Afficher les résultats
cat("Variance empirique totale :", variance_empirique_totale, "\n")
# Variance empirique totale : 0.3815305 
cat("Variance théorique totale :", variance_theorique_totale, "\n")
# Variance théorique totale : 0.3701919 


# Fonctions ---------------------------------------------------------------

## Calcul de la variance théorique de l'effet d'un QTL ----

var_theorique <- function(b, c, q1, q2, phi1){
  
  ## Paramètres :
  p1 <- 1 - q1
  p2 <- 1 - q2
  phi2 <- 1 - phi1
  
  ## Constantes :
  C1 <- (0.5 * p1 * phi1 * q1) + (q1^2 * (phi1 * (1 - phi1)))
  C2 <- (0.5 * p2 * phi2 * q2) + (q2^2 * (phi2 * (1 - phi2)))
  C3 <- q1 * q2 * phi1 * phi2
  
  ## Calcul variance théorique avec les effets B et C : 
  
  var_theo <- (b^2 * C1) + (c^2 * C2) - (2 * b * c * C3)
  
  return(var_theo)
  
}

## Calcul de la variance théorique du score du marqueur (genotype) ----
# codage { 0, 1, 2}
# proba cross 1 : {p1^2, 2p1q1, q1^2}
# proba cross 2 : {p2^2, 2p2q2, q2^2}

var_theo_score <- function(q1, q2, phi1){
  
  ## Paramètres :
  p1 <- 1 - q1
  p2 <- 1 - q2
  phi2 <- 1 - phi1
  
  ## Constantes :
  C1 <- p1 * phi1 + 2 * q1 * (phi1 *(1 - phi1))
  C2 <- p2 * phi2 + 2 * q2 * (phi2 * (1 - phi2))
  C3 <- 8 * q1 * q2 * phi1 * phi2
  
  ## Calcul variance théorique avec les effets B et C : 
  
  var_theo <- 2* q1 * C1 + 2 * q2 * C2 - C3
  
  return(var_theo)
  
}

# test : 
i = 3
essai_theo <- var_theo_score(q1 = tab_param_simu[i,]$q1, 
                             q2 = tab_param_simu[i,]$q2, phi1 = 0.5)

mk_xi <- tab_mk_simu[,i]

geno_bis <- ifelse(mk_xi == "AA", 0, 
                   ifelse((mk_xi == "BB") | (mk_xi == "CC"), 2, 1))

essai_empi <- var(geno_bis)

# Afficher les résultats
cat("Variance empirique  :", essai_empi, "\n")
cat("Variance théorique :", essai_theo, "\n")

# Boucle EFFET pour simuler ---------------------------------------------------

n_mk <- 1000 # nb mk 
N <- 2000 # # Nombre d'individus
# N_cross1 <- 100 # Nombre d'individus cross 1
# N_cross2 <- 100

#  Stocker génotype mk simulé dans une matrice :
tab_mk_simu <- matrix(nrow = N, ncol= n_mk)

# Stocker les paramètres de la simulation :
tab_param_simu <- data.frame(N = rep(N, n_mk),
                             q1 = rep(NA, n_mk),
                             q2 = rep(NA, n_mk), 
                             phi1 = rep(NA, n_mk),
                             b = rep(NA, n_mk),
                             c = rep(NA, n_mk),
                             var_empi = rep(NA, n_mk),
                             var_theo = rep(NA, n_mk))
rownames(tab_param_simu)<-paste0("mk", 1:n_mk)


for (i in 1: n_mk){
  
  print(i)
  
  # Fréquence des allèles :
  ## Tirage dans loi uniforme bornée [0.05 ; 0.5] pour q1 et q2 :
  q1 <- runif(n = 1, min = 0.05, max = 0.5)
  p1 <- 1 - q1
  # q2 <- runif(n = 1, min = 0.05, max = 0.5)
  q2 <- 0
  p2 <- 1 - q2
  
  # Probabilités d'être dans chaque croisement
  phi1 <- 0.5
  phi2 <- 1 - phi1
  
  # Effets pour chaque croisement 
  b <- 3
  c <- b^2
  
  ## après faire une boucle avec 
  # B = 1 / C = 2
  # B = 2 / C = 2
  # B = 3 / C = 3^2
  
  # Simuler les génotypes pour chaque croisement
  genotypes1 <- sample(c("AA", "AB", "BB"), N * phi1, 
                       replace = TRUE, prob = c(p1^2, 2*p1*q1, q1^2))
  genotypes2 <- sample(c("AA", "AC", "CC"), N * phi2, 
                       replace = TRUE, prob = c(p2^2, 2*p2*q2, q2^2))
  genotypes_total <- c(genotypes1, genotypes2)
  # proba sous HW
  
  # Stocker génotype mk simulé dans une matrice :
  tab_mk_simu[,i]<- genotypes_total
  
  # Simuler les phénotypes pour chaque croisement
  phenotypes1 <- ifelse(genotypes1 == "AA", 0, 
                        ifelse(genotypes1 == "BB", b, b/2))
  phenotypes2 <- ifelse(genotypes2 == "AA", 0, 
                        ifelse(genotypes2 == "CC", c, c/2))

  # Calculer la variance globale empirique
  phenotypes_total <- c(phenotypes1, phenotypes2)
  variance_empirique_totale <- var(phenotypes_total)
  
  # Calculer la variance théorique globale
  variance_theorique_totale <- var_theorique(b, c, q1, q2, phi1)
  
  # Stocker paramètres simulation :
  tab_param_simu[i,2:8]<-c(q1, q2, phi1, b, c,
                           variance_empirique_totale,
                           variance_theorique_totale)
  
  # save(tab_mk_simu, file=  "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu.RData")
  # save(tab_param_simu, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu.RData")
}

rownames(tab_mk_simu)<- c(paste0("cross1_", 1:(N*phi1)),
                          paste0("cross2_", 1:(N*phi2)))
colnames(tab_mk_simu)<-paste0("mk", 1:n_mk)

save(tab_mk_simu, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_carre_fixe.RData")
save(tab_param_simu, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_carre_fixe.RData")



# Scatter plot EFFET -----------------------------------------------------------

plot(tab_param_simu$var_empi, tab_param_simu$var_theo,
     pch = 19)
cor(tab_param_simu$var_empi, tab_param_simu$var_theo)

# Les 2 cross segregent : 
## B = 1 // C = 2 
# 0.9541528 pour N = 200 
# 0.9904467 pour N = 1000
# 0.9950611 pour N = 2000

## B = 2 // C = 2 : 
# 0.9333562 pour N = 200
# 0.9839716 pour N = 1000 
# 0.9919822 pour N = 2000

## B = 3 // C = B^2 = 9
# 0.9603404 pour N = 200
# 0.9916013 pour N = 1000
# 0.9957009 pour N = 2000

# En fixant le cross 2 (q2 = 0)
## B = 1 // C = 2 
# 0.9879764 pour N = 200 
# 0.9918092 pour N = 1000
# 0.9960455 pour N = 2000

## B = 2 // C = 2 : 
# 0.9654283 pour N = 200
# 0.9929762 pour N = 1000
# 0.9962355 pour N = 2000 

## B = 3 // C = B^2 = 9
# 0.9655531 pour N = 200
# 0.992406 pour N = 1000
# 0.9963592 pour N = 2000

# Droite de régression :

reg <- lm(var_empi ~ var_theo, data= tab_param_simu)
abline(reg, col="cyan", lty = 2, lwd= 2)


# Simuler variance --------------------------------------------------------

load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc.RData")
# c'est des simulations avec B = C 
# même si les autres param sont pas influencés pour ca 
# c'est juste histoire de...

i = 3

## 1) Fixer sigma^2G = 50 ----

sigma_G <- 50

## 2) Sélectionner 1 marqueur xi ----

mk_xi <- tab_mk_simu[,i]


## 3) Tirer effet du marqueur ----

# On pose b = c
# D'où, sigma^2G = b^2*C1 + c^2 * C2 - 2 * b * c * C3
#                = b^2*C1 + b^2 * C2 - 2 * b * b * C3
#                = b^2 (C1 + C2 - 2 * C3)
# On a, b = +/- racine(sigma^2G / (C1+C2-2*C3))
# Donc beta = racine(50 / (C1+C2-2*C3))

# Paramètres pour calculer les constantes :
q1_bis <- tab_param_simu[i,]$q1
p1_bis <- 1 - q1_bis
q2_bis <- tab_param_simu[i,]$q2
p2_bis <- 1 - q2_bis
phi1_bis <- tab_param_simu[i,]$phi1
phi2_bis <- 1 - phi1_bis

C1_bis <- (0.5 * p1_bis * phi1_bis * q1_bis) + (q1_bis^2 * (phi1_bis * (1 - phi1_bis)))
C2_bis <- (0.5 * p2_bis * phi2_bis * q2_bis) + (q2_bis^2 * (phi2_bis * (1 - phi2_bis)))
C3_bis <- q1_bis * q2_bis * phi1_bis * phi2_bis

beta_i <- sqrt(sigma_G/(C1_bis+ C2_bis - 2*C3_bis))

## 4) Calculer phénotype à ce marqueur : xi * beta ----
 
# recodage en 0, 1, 2 
geno_bis <- ifelse(mk_xi == "AA", 0, 
                   ifelse((mk_xi == "BB") | (mk_xi == "CC"), 2, 1))

pheno_bis <- geno_bis * beta_i

## 5) Calculer la variance du phénotype simulé ----

var(pheno_bis)
# 184.2707
# 191.5543
# 219.9767

# Est-ce que cette variance est égale à sigma^2G = 50 ?
# ... pas du tout ...

## Boucle PARAM SIMU SIGMA^2G = 50 : ----

# N = 200
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc.RData")

# N = 1000
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N1000_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N1000_bc.RData")

# N = 2000
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_bc.RData")

# N = 200 FIXE
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc_fixe.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc_fixe.RData")

# N = 1000 FIXE
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N1000_bc_fixe.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N1000_bc_fixe.RData")

# N = 2000 FIXE
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_bc_fixe.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_bc_fixe.RData")



sigma_G <- 50

tab_res_variances <- data.frame(var_sim = rep(NA,nrow(tab_param_simu)))

for (i in 1:nrow(tab_param_simu)){
  
  ## Sélectionner 1 marqueur xi
  mk_xi <- tab_mk_simu[,i]
  
  ## Tirer effet du marqueur
  
  # Paramètres pour calculer les constantes :
  q1_bis <- tab_param_simu[i,]$q1
  p1_bis <- 1 - q1_bis
  q2_bis <- tab_param_simu[i,]$q2
  p2_bis <- 1 - q2_bis
  phi1_bis <- tab_param_simu[i,]$phi1
  phi2_bis <- 1 - phi1_bis
  
  C1_bis <- (0.5 * p1_bis * phi1_bis * q1_bis) + (q1_bis^2 * (phi1_bis * (1 - phi1_bis)))
  C2_bis <- (0.5 * p2_bis * phi2_bis * q2_bis) + (q2_bis^2 * (phi2_bis * (1 - phi2_bis)))
  C3_bis <- q1_bis * q2_bis * phi1_bis * phi2_bis
  
  # u <- sample(c(1, -1), size = 1)
  # c'est un calcul de variance donc blc de si c'est positif ou neg
  
  beta_i <- sqrt(sigma_G/(C1_bis+ C2_bis - 2*C3_bis))
  
  ## Calculer phénotype à ce marqueur
  
  geno_bis <- ifelse(mk_xi == "AA", 0, 
                     ifelse((mk_xi == "BB") | (mk_xi == "CC"), 2, 1))
  
  
  # Standardiser le génotype
  geno_bis <- (geno_bis - mean(geno_bis)) / sd(geno_bis)
  
  pheno_bis <- geno_bis * beta_i
  # pheno_bis <- 0.5 * (geno_bis * beta_i)
  
  ## Calculer la variance du phénotype simulé
  
  tab_res_variances[i,1] <- var(pheno_bis)
  
  
  
}


# save(tab_res_variances, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_res_variances_1000_N2000_bc_fixe.RData")

boxplot(tab_res_variances$var_sim)
summary(tab_res_variances$var_sim)

# Pour N = 200 en simulant :

## SANS 1/2 : c'est grave surestimé 
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 141.5   187.8   198.9   198.9   209.7   270.1 

## AVEC 1/2 : What the hell, ça marche !?!?
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 35.38   46.96   49.72   49.71   52.41   67.52 

# Pour N = 1000 en simulant : (tjrs avec 1/2)

#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 41.58   48.67   50.03   50.04   51.38   61.54 

# Pour N = 2000 en simulant : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.28   49.07   50.01   50.05   51.01   56.11 

## Avec 1 cross fixé : 

# N = 200 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.29   44.95   49.49   49.70   54.11   80.15 



# BOUCLE Variance théorique du score des marqueurs -------------------------------

# load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_essai_fixe.RData")
# load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_essai_fixe.RData")

# load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N1000_essai_fixe.RData")
# load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N1000_essai_fixe.RData")

load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_essai_fixe.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_essai_fixe.RData")

n_mk <- 1000 # nb mk 
N <- 2000 # Nombre d'individus

# Stocker les paramètres de la simulation :
tab_param_score <- data.frame(N = rep(N, n_mk),
                             q1 = rep(NA, n_mk),
                             q2 = rep(NA, n_mk), 
                             phi1 = rep(NA, n_mk),
                             var_empi = rep(NA, n_mk),
                             var_theo = rep(NA, n_mk))
rownames(tab_param_score)<-paste0("mk", 1:n_mk)


for (i in 1: n_mk){
  
  print(i)
  
  # Fréquence des allèles :
  q1 <- tab_param_simu[i,]$q1
  p1 <- 1 - q1
  q2 = tab_param_simu[i,]$q2
  p2 <- 1 - q2
  
  # Probabilités d'être dans chaque croisement
  phi1 <- 0.5
  phi2 <- 1 - phi1

  # Récuperer génotypes des mk simulés avant : 
  
  mk_xi <- tab_mk_simu[,i]
  
  geno_bis <- ifelse(mk_xi == "AA", 0, 
                     ifelse((mk_xi == "BB") | (mk_xi == "CC"), 2, 1))
  
  
  # Calculer la variance globale empirique
  variance_empirique_totale <- var(geno_bis)
  
  # Calculer la variance théorique globale
  variance_theorique_totale <- var_theo_score(q1, 
                               q2 , phi1)

  
  # Stocker paramètres simulation :
  tab_param_score [i,2:6]<-c(q1, q2, phi1,
                           variance_empirique_totale,
                           variance_theorique_totale)
  
  # save(tab_mk_simu, file=  "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu.RData")
  # save(tab_param_simu, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu.RData")
}


save(tab_param_score, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_score_1000_N2000_fixe.RData")


# Scatter plot SCORE -----------------------------------------------------------

plot(tab_param_score$var_empi, tab_param_score$var_theo,
     pch = 19)
cor(tab_param_score$var_empi, tab_param_score$var_theo)


# Les 2 cross segregent : 
# 0.9313579 pour N = 200 
# 0.9849681 pour N = 1000
# 0.9921671 pour N = 2000


# En fixant le cross 2 (q2 = 0)
# 0.962897 pour N = 200 
# 0.9918092 pour N = 1000
# 0.9960455 pour N = 2000


# Droite de régression :

reg <- lm(var_empi ~ var_theo, data= tab_param_score)
abline(reg, col="cyan", lty = 2, lwd= 2)


# c = sb : Variance effet QTL : c = sb ---------------------------------------------

# Pop F2 avec fq all équilibrées

library(simcross)
library(mppSim)
library(mppR)

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

# BOUCLE simulation of the variance of a single QTL ----
X1 <- geno[1:100, ]
X2 <- geno[101:200, ]
X <- rbind(X1, X2)
# Xs <- scale_center_X(X = X)

# var_s <- var_ns <- rep(NA, 100)
var_s <- rep(NA, 1000)

# set the genetic variance
sigma_G <- 50

for(i in 1:1000){
  
  # pos_i <- sample(1:ncol(X), size = 1)
  pos_i <- i
  
  X1_sel <- X1[, pos_i, drop = FALSE]
  X2_sel <- X2[, pos_i, drop = FALSE]
  X_sel <- rbind(X1_sel, X2_sel)
  # X_sel_s <- Xs[, pos_i, drop = FALSE]
  
  # calculate the marker frequency
  q1 <- apply(X1_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x)))
  q2 <- apply(X2_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x)))
  
  p1 <- 1 - q1
  p2 <- 1 - q2
  
  # Probabilités d'être dans chaque croisement
  phi1 <- 0.5
  phi2 <- 1 - phi1
  
  # Constantes
  C1 <- (0.5 * p1 * phi1 * q1) + (q1^2 * (phi1 * (1 - phi1)))
  C2 <- (0.5 * p2 * phi2 * q2) + (q2^2 * (phi2 * (1 - phi2)))
  C12 <- 2 * q1 * q2 * phi1 * phi2
  
  # s <- 0.6^2
  
  # determine the QTL effect
  
  f <- 1/2
  
  u <- sample(c(1, -1), size = 1)
  
  beta1_i <- u * f *  sqrt(sigma_G/(C1 + s^2 * C2 - s* C12))
  # effet cross 1

  beta2_i <- s * beta1_i
  # effet cross 2
  
  # beta_i <- u * f *  sqrt(sigma_G/(C1 + C2 - C12))
  
  pheno1 <- X1_sel * beta1_i
  pheno2 <- X2_sel * beta2_i
  
  pheno <- rbind(pheno1,pheno2)
  
  # pheno <- X_sel * beta_i
  
  var_s[i] <- var(pheno)
  
}




res <- data.frame(var_s)
res_effets <- res

save(res_effets, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/res_effets.RData")

colMeans(res)


boxplot(res)

abline(h = 50, col = "green")
summary(res)
# abline(h = 100, col = "blue")

# Avec s = 0.9^2
# Avec s = 0.6^2 ça marche aussi 


# c = sb : FQ ALL DESEQ ---------------------------------------------------

# N = 200
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc.RData")

# N = 2000
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_bc.RData")

# FIXE : 
# N = 200
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc_fixe.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc_fixe.RData")

# N = 2000
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_bc.RData")


tab_mk_cross1 <- tab_mk_simu[1:(0.5 * nrow(tab_mk_simu)), ]
tab_mk_cross2 <- tab_mk_simu[(0.5 * nrow(tab_mk_simu)+1):nrow(tab_mk_simu), ]


sigma_G <- 50

tab_res_variances <- data.frame(var_sim = rep(NA,nrow(tab_param_simu)))

for (i in 1:nrow(tab_param_simu)){
  
  ## Sélectionner 1 marqueur xi
  mk_xi_cross1 <- matrix(tab_mk_cross1[,i])
  mk_xi_cross2 <- matrix(tab_mk_cross2[,i])
  
  # mk_xi <- rbind(mk_xi_cross1, mk_xi_cross2)
  # rownames(mk_xi) <-rownames(tab_mk_simu)
  
  ## Tirer effet du marqueur
  
  # Paramètres pour calculer les constantes :
  q1 <- tab_param_simu[i,]$q1
  p1 <- 1 - q1
  q2 <- tab_param_simu[i,]$q2
  p2 <- 1 - q2
  phi1 <- tab_param_simu[i,]$phi1
  phi2 <- 1 - phi1
  
  C1 <- (0.5 * p1 * phi1 * q1) + (q1^2 * (phi1 * (1 - phi1)))
  C2 <- (0.5 * p2 * phi2 * q2) + (q2^2 * (phi2 * (1 - phi2)))
  C3 <- 2 * q1 * q2 * phi1 * phi2
  
  s <- 0.9^2
  
  u <- sample(c(1, -1), size = 1)
  # c'est un calcul de variance donc blc de si c'est positif ou neg
  
  f <- 1/2
  
  # beta_i <- sqrt(sigma_G/(C1_bis+ C2_bis - C3_bis))
  
  beta1_i <- u * f *  sqrt(sigma_G/(C1 + s^2 * C2 - s* C3))
  # effet cross 1
  
  beta2_i <- s * beta1_i
  
  ## Calculer phénotype à ce marqueur
  
  geno1 <- ifelse(mk_xi_cross1 == "AA", 0, 
                     ifelse((mk_xi_cross1 == "BB") | (mk_xi_cross1 == "CC"), 2, 1))
  
  geno2 <- ifelse(mk_xi_cross2 == "AA", 0, 
                  ifelse((mk_xi_cross2 == "BB") | (mk_xi_cross2 == "CC"), 2, 1))
  
  
  pheno1 <- geno1 * beta1_i
  pheno2 <- geno2 * beta2_i
  
  pheno <- rbind(pheno1, pheno2)
  # pheno_bis <- 0.5 * (geno_bis * beta_i)
  
  ## Calculer la variance du phénotype simulé
  
  tab_res_variances[i,1] <- var(pheno)
  
  
  
}

save(tab_res_variances, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_res_var_1000_N200_fixe.RData")

boxplot(tab_res_variances$var_sim)
abline(h=50, col= 'cyan', lwd = 2, lty= 2)

summary(tab_res_variances$var_sim)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.53   46.73   49.64   49.73   52.43   68.66 



# Plusieurs marqueurs  ----------------------------------------------------

# Avec 10 marqueurs ayant des proportions de variance expliquée
# par le QTL décroissante :

# Genre : ki = 1 / i^2 
# puis normaliser pour que la somme = 1 

## B = C ----

# N = 200
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc.RData")

# N = 2000
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_bc.RData")


sigma_G <- 50
nb_mk <- 1000

# Définir les proportions décroissantes k_i
k <- 1 / (1:nb_mk)^2
k <- k / sum(k) # Normaliser pour que la somme soit 1


# k <- 0.5 * 0.5^(1:nb_mk)
# k <- k / sum(k) 
# 
# k<- 1 - 0.1*(1:nb_mk)
# k <- k / sum(k) 

n_rep <- 1

tab_somme <-data.frame(somme = rep(NA, n_rep))

for (j in 1:n_rep){
  
  tab_res_variances <- data.frame(var_sim = rep(NA,nrow(tab_param_simu)))
  
  for (i in 1:nb_mk){
    
    ## Sélectionner 1 marqueur xi
    mk_xi <- tab_mk_simu[,i]
    # nom_mk_xi <- sample(colnames(tab_mk_simu), 1)
    # mk_xi<- tab_mk_simu[,nom_mk_xi]
    
    ## Tirer effet du marqueur
    
    # Paramètres pour calculer les constantes :
    q1_bis <- tab_param_simu[i,]$q1
    p1_bis <- 1 - q1_bis
    q2_bis <- tab_param_simu[i,]$q2
    p2_bis <- 1 - q2_bis
    phi1_bis <- tab_param_simu[i,]$phi1
    phi2_bis <- 1 - phi1_bis
    
    C1_bis <- (0.5 * p1_bis * phi1_bis * q1_bis) + (q1_bis^2 * (phi1_bis * (1 - phi1_bis)))
    C2_bis <- (0.5 * p2_bis * phi2_bis * q2_bis) + (q2_bis^2 * (phi2_bis * (1 - phi2_bis)))
    C3_bis <- 2 * q1_bis * q2_bis * phi1_bis * phi2_bis
    
    # u <- sample(c(1, -1), size = 1)
    # c'est un calcul de variance donc blc de si c'est positif ou neg
    
    beta_i <- 0.5 * sqrt((k[i] * sigma_G)/(C1_bis+ C2_bis - C3_bis))
    
    ## Calculer phénotype à ce marqueur
    
    geno_bis <- ifelse(mk_xi == "AA", 0, 
                       ifelse((mk_xi == "BB") | (mk_xi == "CC"), 2, 1))
    
    pheno_bis <- geno_bis * beta_i
    
    ## Calculer la variance du phénotype simulé
    
    tab_res_variances[i,1] <- var(pheno_bis)
    
    
    
  }
  
  tab_somme[j,1] <- sum(tab_res_variances$var_sim, na.rm = TRUE)
  
}



# La somme des 10 mk doit valoir sigma^2_G = 50 : 

save(tab_res_variances, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_somme_1000_N2000.RData")

boxplot(tab_somme$somme)
abline(h= 50, col = "cyan", lty = 2)
summary(tab_somme$somme)


## POP avec N=200 : 
# Pour 1000 rep et ki = 1/i^2: 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 24.04   42.48   48.67   47.22   52.93   64.41 

# Pour 1000 rep et ki = 0.5 * 0.5^i
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 24.15   42.28   47.36   46.70   51.90   64.39

# Pour 2000 rep : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 22.54   43.17   47.89   46.98   51.63   65.04

## POP avec N=2000

# Pour 1000 rep et ki = 1/i^2: 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 26.95   52.95   59.02   57.76   63.78   71.68 

# 2000 rep :
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 30.56   53.22   59.33   58.04   64.33   71.80 

# 2000 Rep l'autre ki : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 30.27   56.73   62.92   61.67   67.92   76.85 

# N = 200 / 1000 mk, 1000 rep :
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 20.93   40.09   45.63   44.47   49.89   61.54 

# N = 2000 / 1000 mk, 1000rep
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 27.64   49.70   55.89   54.56   60.35   68.48 

## C = SB ----

# N = 200
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N200_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N200_bc.RData")

# N = 2000
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_mk_simu_1000_N2000_bc.RData")
load("~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_param_simu_1000_N2000_bc.RData")


tab_mk_cross1 <- tab_mk_simu[1:(0.5 * nrow(tab_mk_simu)), ]
tab_mk_cross2 <- tab_mk_simu[(0.5 * nrow(tab_mk_simu)+1):nrow(tab_mk_simu), ]



# 3 CROSS ----

## Pop F2 avec fq all équilibrées ----

library(simcross)
library(mppSim)
library(mppR)
library(LDcorSV)
library(ggplot2)

# sub-functions

scale_center_X <- function(X){
  
  p <- apply(X, 2, FUN = function(x) sum(x)/(2*length(x)))
  P <- rep(1, nrow(X)) %*% (2*t(p))
  Z <- X - P
  den <- sqrt(2 * p * (1 - p))
  Z <- t(t(Z) / den)
  # equivalent to Z %*% diag(1 / den)
  
  return(Z)
  
}

# simulation of a 8 crosses NAM population (9 parents)

data("geno_par")
rownames(geno_par) <- paste0('P', 1:nrow(geno_par))

data("EUNAM_map")
rownames(map) <- map[, 1]

# NAM crossing scheme with 9 parents
cross_scheme <- cross_scheme_NAM(9)

geno <- sim_mpp_cross(geno_par = geno_par, map = map,
                      cross_scheme = cross_scheme,
                      geno_score_format = "numeric_IBD")

### BOUCLE simulation of the variance of a single QTL ----
X1 <- geno[1:100, ]
X2 <- geno[101:200, ]
X3 <- geno[201:300, ]
X <- rbind(X1, X2, X3)
# Xs <- scale_center_X(X = X)

# var_s <- var_ns <- rep(NA, 100)
var_s <- rep(NA, 1000)

# set the genetic variance
sigma_G <- 50

for(i in 1:1000){
  
  # pos_i <- sample(1:ncol(X), size = 1)
  pos_i <- i
  
  X1_sel <- X1[, pos_i, drop = FALSE]
  X2_sel <- X2[, pos_i, drop = FALSE]
  X3_sel <- X3[, pos_i, drop = FALSE]
  X_sel <- rbind(X1_sel, X2_sel, X3_sel)
  # X_sel_s <- Xs[, pos_i, drop = FALSE]
  
  # calculate the marker frequency
  q1 <- as.numeric(apply(X1_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x))))
  q2 <- as.numeric(apply(X2_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x))))
  q3 <- as.numeric(apply(X3_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x))))
  
  
  p1 <- 1 - q1
  p2 <- 1 - q2
  p3 <- 1 - q3
  
  # Probabilités d'être dans chaque croisement
  phi1 <- 1/3
  phi2 <- 1/3
  phi3 <- 1/3
  
  # Constantes
  C1 <- (0.5 * p1 * phi1 * q1) + (q1^2 * (phi1 * (1 - phi1)))
  C2 <- (0.5 * p2 * phi2 * q2) + (q2^2 * (phi2 * (1 - phi2)))
  C3 <- (0.5 * p3 * phi3 * q3) + (q3^2 * (phi3 * (1 - phi3)))
  C12 <- 2 * q1 * q2 * phi1 * phi2
  C13 <- 2 * q1 * q3 * phi1 * phi3
  C23 <- 2 * q2 * q3 * phi2 * phi3
  
  # s <- 0.9^2
  
  # determine the QTL effect
  
  f <- 1/2
  
  u <- sample(c(1, -1), size = 1)
  
  # beta_i <- u * f * sqrt((sigma_G)/((C1 + C2 + C3) -
  # ( C12 + C13 + C23)))

  beta1_i <- u * f *  sqrt(sigma_G/( (C1 + s^2 * C2 + s^4* C3) -
                                     (s * C12 + s^2 * C13 + s^3 * C23)))
  # # effet cross 1
  #
  beta2_i <- s * beta1_i
  # # effet cross 2
  #
  beta3_i <- s^2 * beta1_i
  #
  pheno1 <- X1_sel * beta1_i
  pheno2 <- X2_sel * beta2_i
  pheno3 <- X3_sel * beta3_i

  pheno <- rbind(pheno1,pheno2, pheno3)

  # pheno <- X_sel * beta_i
  
  var_s[i] <- var(pheno)
  
}
res <- data.frame(var_s)

res_3cross_egal<-res

# save(res_3cross_sbcd, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/res_3cross_sbcd.RData")
save(res_3cross_egal, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/res_3cross_egal.RData")

boxplot(res)
abline(h = 50, col = "cyan")
summary(res$var_s)

## pour b = c = d
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 42.88   48.33   50.44   50.41   52.27   56.38 

# pour effet lié : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 41.85   48.41   50.38   50.48   52.39   55.96 

## FQ DESEQ ----

### simuler fq all deseq pour 3 cross
### verifier b = c =d 
### et sbcd



tab_mk_cross1 <- tab_mk_simu[1:(0.5 * nrow(tab_mk_simu)), ]
tab_mk_cross2 <- tab_mk_simu[(0.5 * nrow(tab_mk_simu)+1):nrow(tab_mk_simu), ]
tab_mk_cross3 <- tab_mk_simu[(0.5 * nrow(tab_mk_simu)+1):nrow(tab_mk_simu), ]


sigma_G <- 50

tab_res_variances <- data.frame(var_sim = rep(NA,nrow(tab_param_simu)))

for (i in 1:nrow(tab_param_simu)){
  
  ## Sélectionner 1 marqueur xi
  mk_xi_cross1 <- matrix(tab_mk_cross1[,i])
  mk_xi_cross2 <- matrix(tab_mk_cross2[,i])
  
  # mk_xi <- rbind(mk_xi_cross1, mk_xi_cross2)
  # rownames(mk_xi) <-rownames(tab_mk_simu)
  
  ## Tirer effet du marqueur
  
  # Paramètres pour calculer les constantes :
  q1 <- tab_param_simu[i,]$q1
  p1 <- 1 - q1
  q2 <- tab_param_simu[i,]$q2
  p2 <- 1 - q2
  phi1 <- tab_param_simu[i,]$phi1
  phi2 <- 1 - phi1
  
  C1 <- (0.5 * p1 * phi1 * q1) + (q1^2 * (phi1 * (1 - phi1)))
  C2 <- (0.5 * p2 * phi2 * q2) + (q2^2 * (phi2 * (1 - phi2)))
  C3 <- 2 * q1 * q2 * phi1 * phi2
  
  s <- 0.9^2
  
  u <- sample(c(1, -1), size = 1)
  # c'est un calcul de variance donc blc de si c'est positif ou neg
  
  f <- 1/2
  
  # beta_i <- sqrt(sigma_G/(C1_bis+ C2_bis - C3_bis))
  
  beta1_i <- u * f *  sqrt(sigma_G/(C1 + s^2 * C2 - s* C3))
  # effet cross 1
  
  beta2_i <- s * beta1_i
  
  ## Calculer phénotype à ce marqueur
  
  geno1 <- ifelse(mk_xi_cross1 == "AA", 0, 
                  ifelse((mk_xi_cross1 == "BB") | (mk_xi_cross1 == "CC"), 2, 1))
  
  geno2 <- ifelse(mk_xi_cross2 == "AA", 0, 
                  ifelse((mk_xi_cross2 == "BB") | (mk_xi_cross2 == "CC"), 2, 1))
  
  
  pheno1 <- geno1 * beta1_i
  pheno2 <- geno2 * beta2_i
  
  pheno <- rbind(pheno1, pheno2)
  # pheno_bis <- 0.5 * (geno_bis * beta_i)
  
  ## Calculer la variance du phénotype simulé
  
  tab_res_variances[i,1] <- var(pheno)
  
  
  
}

save(tab_res_variances, file = "~/M2 DATA ANALYST/STAGE/BCNAM_Data/Scripts_Data/Simulations/tab_res_var_1000_N200_fixe.RData")

boxplot(tab_res_variances$var_sim)
abline(h=50, col= 'cyan', lwd = 2, lty= 2)

summary(tab_res_variances$var_sim)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.53   46.73   49.64   49.73   52.43   68.66 




#### GENERALISATION : 




# GENERALISATION ----------------------------------------------------------

## Fq équilibrées : ----

### effets égaux ----

#### fonctions ----
# Constantes interactions : 

mc_const_inter <- function(ind_i, ind_j){
  
  C_i_j <- 2 * ind_i$q_i * ind_j$q_i * ind_i$phi_i * ind_j$phi_i
  return(C_i_j)
}

# Dénominateur effets égaux (utilise const_inter)
mc_calcul_denom <- function(nc, tab_param_cross){
  
  # partie 1 :
  somme_p1 <- tab_param_cross[[1]]$C_i
  for (i in 2:nc){
    # print(paste("somme = ", somme_p1))
    somme_p1 <- somme_p1 + tab_param_cross[[i]]$C_i
    # print(paste("somme = ", somme_p1))
  }
  
  # partie 2 : 
  somme_p2 <- 0
  for (i in 1:nc){
    # print (paste("i =", i))
    for (j in (i+1):nc){
      # print(paste("j=", j))
      if ( (i!= j) & (j<=length(tab_param_cross))){
        # print("coucou")
        # print(paste("somme = ", somme_p2))
        inter <- const_inter(tab_param_cross[[i]], tab_param_cross[[j]])
        somme_p2 <- somme_p2 + inter
        # print(paste("somme = ", somme_p2))
      }
    }
  }
  
  denom <- somme_p1 - somme_p2
  return(denom)
  
}

# Récupérer paramètres croisements :
mc_param_cross <- function (geno, nc, pos_i){
  
  cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)
  cr_id <-unique(cr_ind)
  
  tab_param_cross <- vector(mode = 'list', length = nc)
  
  for (i in 1:nc){
    
    geno_i<-geno[cr_ind %in% cr_id[i], ]
    
    geno_i_sel <- geno_i[, pos_i, drop = FALSE]
    q_i <- as.numeric(apply(geno_i_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x))))
    
    p_i <- 1 - q_i
    
    phi_i <- (nrow(geno)/nc) / nrow(geno)
    
    # Constante Ci : 
    
    C_i <- (0.5 * p_i * q_i * phi_i) + (q_i^2 * (phi_i * (1 - phi_i)))
    
    list_i <- list(geno_i = geno_i, 
                   geno_i_sel = geno_i_sel, 
                   q_i = q_i,
                   p_i = p_i,
                   phi_i = phi_i, 
                   C_i = C_i)
    
    
    tab_param_cross[[i]]<-list_i
    
    
  }
  
  return(tab_param_cross)
  
}




#### Boucle ---- 

n_rep <-100 # nb de rep, 1 mk à chaque fois

n_cross <-2

var_s <- rep(NA, n_rep)

sigma_G <- 50

for (k in 1:n_rep){
  
  
  pos_i <- sample(1:ncol(geno), size = 1)
  # pos_i <- k
  

  tab_param_cross <- mc_param_cross(geno, nc = n_cross, pos_i )
  
  denom <- mc_calcul_denom(nc = n_cross, 
                        tab_param_cross = tab_param_cross)
  
  
  # determine the QTL effect

  f <- 1/2

  u <- sample(c(1, -1), size = 1)

  beta_i <- u * f * sqrt((sigma_G)/ denom )


  geno_sel<-data.frame(tab_param_cross[[1]]$geno_i_sel)
  for (g in 2:length(tab_param_cross)){
    geno_sel <- rbind(geno_sel, tab_param_cross[[g]]$geno_i_sel)
  }

  pheno <- geno_sel * beta_i
  
  var_s[k] <- var(pheno)
  
  
}

var_res<-data.frame(var_s)
boxplot(var_res)
abline(h=50, col= "cyan", lty= 2)
summary(var_res$var_s)

# 100 rep : 2 cross
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 43.85   49.72   51.82   51.67   55.08   57.34 

# 100 rep : 3 cross 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.97   50.53   52.07   51.91   53.42   58.12 

# 100 rep : 4 cross
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.27   49.68   51.29   51.11   52.49   57.13

# 100 rep : 5 cross
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.92   49.23   51.17   50.99   52.89   54.59 

# 100 rep : 6 cross
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.51   48.82   50.47   50.33   51.86   54.07 

# 100 rep : 7 cross 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.98   49.07   50.47   50.27   51.67   54.77 

# 100 rep : 8 cross
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.98   49.55   50.55   50.42   51.35   54.48






# est-ce que ça va vrmt se rapprocher de 50 + on augmente nb
# cross (ou nb d'individus ? ), ou c'est genre avec 
# 700 - 800 indiv 
# que c'est le mieux


# Je sais pas si c'est "mieux" parce que y'a + d'individus 
# comme + de croisement (et que c'est effet constant)
# ou bien c'est vrmt en ajoutant des croisements ?
# mais du coup, t'façon là c'est fq alléliques équilibrées
# faut que je fasse la fonction pour effet liés 
# pour voir si c'est qqchose de similaire ou pas ?


### effets liés ----

#### fonctions ----

# Constantes interactions : 

mc_const_inter <- function(ind_i, ind_j){
  
  C_i_j <- 2 * ind_i$q_i * ind_j$q_i * ind_i$phi_i * ind_j$phi_i
  return(C_i_j)
}

# Dénominateur effets égaux (utilise const_inter)
mc_calcul_denom_s <- function(nc, tab_param_cross, s){
  
  # partie 1 :
  somme_p1 <- 0
  for (i in 1:nc){
    # print(paste("somme = ", somme_p1))
    
    somme_p1 <- somme_p1 + s^(2*(i-1)) *  tab_param_cross[[i]]$C_i
    
    # print(paste("somme = ", somme_p1))
  }
  
  # partie 2 : 
  somme_p2 <- 0
  for (i in 1:nc){
    # print (paste("i =", i))
    for (j in (i+1):nc){
      # print(paste("j=", j))
      if ( (i!= j) & (j<=length(tab_param_cross))){
        # print("coucou")
        # print(paste("somme = ", somme_p2))
        inter <- const_inter(tab_param_cross[[i]], tab_param_cross[[j]])
        somme_p2 <- somme_p2 + s^(i+j-2) * inter
        # print(paste("somme = ", somme_p2))
      }
    }
  }
  
  denom <- somme_p1 - somme_p2
  return(denom)
  
}

# Récupérer paramètres croisements :
mc_param_cross <- function (geno, nc, pos_i){
  
  cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)
  cr_id <-unique(cr_ind)
  
  # nc <- length(cr_id)
  
  tab_param_cross <- vector(mode = 'list', length = nc)
  
  for (i in 1:nc){
    
    geno_i<-geno[cr_ind %in% cr_id[i], ]
    
    geno_i_sel <- geno_i[, pos_i, drop = FALSE]
    q_i <- as.numeric(apply(geno_i_sel, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x))))
    
    p_i <- 1 - q_i
    
    phi_i <- (nrow(geno)/nc) / nrow(geno)
    
    # Constante Ci : 
    
    C_i <- (0.5 * p_i * q_i * phi_i) + (q_i^2 * (phi_i * (1 - phi_i)))
    
    list_i <- list(geno_i = geno_i, 
                   geno_i_sel = geno_i_sel, 
                   q_i = q_i,
                   p_i = p_i,
                   phi_i = phi_i, 
                   C_i = C_i)
    
    
    tab_param_cross[[i]]<-list_i
    
    
  }
  
  return(tab_param_cross)
  
}



#### Boucle ----

n_rep <-100 # nb de rep, 1 mk à chaque fois

n_cross<- 2


# cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)
# cr_id <-unique(cr_ind)

var_s <- rep(NA, n_rep)

sigma_G <- 50

s <- 0.9^2

for (k in 1:n_rep){
  
  
  pos_i <- sample(1:ncol(geno), size = 1)
  # pos_i <- k
  
  tab_param_cross <- mc_param_cross(geno, nc = n_cross,
                                    pos_i )
  
  denom <- mc_calcul_denom_s(nc = n_cross, 
                        tab_param_cross = tab_param_cross, s = s)
  
  # determine the QTL effect
  
  f <- 1/2
  
  u <- sample(c(1, -1), size = 1)
  
  beta1_i <- u * f * sqrt((sigma_G)/ denom)
  pheno_1 <- tab_param_cross[[1]]$geno_i_sel * beta1_i 
  
  mespheno <- vector(mode='list', length = n_cross)
  mespheno[[1]] <- pheno_1
  
  ## pour les autres betas apres 
  for (r in 2:n_cross){
    
    facteur_s <-s^(r-1)
    
    beta_i<- beta1_i * facteur_s
    pheno_i <- tab_param_cross[[r]]$geno_i_sel * beta_i
    mespheno[[r]] <- pheno_i
  }

  phenob <-mespheno[[1]]
  for (r in 2:n_cross){
    phenob <- rbind(phenob, mespheno[[r]])
  }
  
  
  var_s[k] <- var(phenob)
  
}

var_res<-data.frame(var_s)
boxplot(var_res)
abline(h=50, col= "cyan", lty= 2)
summary(var_res$var_s)

## 100 rep : 2 cross
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 41.72   48.05   51.88   51.07   53.73   57.37

# 3 cross : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.81   49.98   51.98   51.40   53.15   55.77 

# 4 cross :
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.57   49.57   50.64   50.75   52.23   55.04 

# 5 cross : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.66   49.94   51.30   51.11   52.56   54.64 

# 6 cross : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.18   49.23   50.92   50.52   52.04   53.93 

# 7 cross :
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.42   49.90   51.22   50.74   51.95   54.53 

# 8 cross : 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.29   49.38   50.74   50.35   51.63   53.17 


# PLUSIEURS MK  -----------------------------------------------------------

## EGAUX ----

###  fonctions ----
mc_allele_effect_equal <- function(k, nc, sigma_G, denom){
  
  f <- 1/2
  u <- sample(c(1, -1), size = 1)
  
  # l'effet égal pour tous les cross :
  beta_i <- u * f * sqrt((k* sigma_G)/ denom)
  
  beta_effets<- vector(mode = "list", length = nc)
  for (i in 1:nc){
    beta_effets[[i]]<- beta_i
  }
  return(beta_effets)
}

mc_calcul_pheno <- function(betas, nc, tab_param_cross){
  
  phenos <- vector(mode='list', length = nc)
  
  for (i in 1:nc){
    
    pheno_i <- tab_param_cross[[i]]$geno_i_sel * betas[[i]]
    phenos[[i]] <- pheno_i
  }
  
  phenob <-phenos[[1]]
  for (i in 2:nc){
    phenob <- rbind(phenob, phenos[[i]])
  }
  
  return(phenob)
}

### BOUCLE -----

n_rep <-100 # nb de rep
n_cross <-2
nb_mk <- 10

# Définir les proportions décroissantes k_i
k_i <- 1 / (1:nb_mk)^2
k_i <- k_i / sum(k_i) # Normaliser pour que la somme soit 1

sigma_G <- 50


tab_somme <-data.frame(somme = rep(NA, n_rep))

for (o in 1:n_rep){
  print(o)
  
  var_s <- rep(NA, n_rep)
  
  for (k in 1:nb_mk){
    
    
    pos_i <- sample(1:ncol(geno), size = 1)
    # pos_i <- k
    
    tab_param_cross <- mc_param_cross(geno, nc = n_cross,pos_i )
    
    denom <- mc_calcul_denom(nc = n_cross, 
                          tab_param_cross = tab_param_cross)
    
    
    # determine the QTL effect
    beta <- mc_allele_effect_equal(k = k_i[k], nc = n_cross, 
                                     sigma_G = sigma_G, 
                                     denom = denom)
    # calcule des phenotypes
    pheno_eq <- mc_calcul_pheno(betas = beta, 
                                     nc = n_cross, 
                                     tab_param_cross = tab_param_cross)
    
    
    var_s[k] <- var(pheno_eq)
    
    var_res<-data.frame(var_s)
    
    
  }
  
  tab_somme[o,1] <- sum(var_res$var_s, na.rm = TRUE)
  
  
  }



# tab_somme

## 2 CROSS : 

# pour 1000 rep avec 10 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 43.05   50.00   51.98   51.65   53.52   57.05 

# pour 500 rep avec 100 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.   
# 45.82   49.75   51.72   51.52   53.27   56.98     

# pour 100 Rep avec 500 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.88   50.56   52.21   51.94   53.67   56.13 


## 5 CROSS :

# pour 100 rep avec 10 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.19   50.06   51.06   50.92   51.96   54.26 

# pour 100 rep avec 100 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 48.23   50.63   51.61   51.41   52.29   53.65 

# pour 100 rep avec 500 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.96   49.95   51.23   51.21   52.51   53.99 

## 8 CROSS: 

# pour 100 rep avec 10 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.09   49.39   50.37   50.32   51.06   52.79 

# pour 100 rep avec 100 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.60   49.52   50.50   50.35   51.01   52.54 

# pour 100 rep avec 500 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.94   49.89   50.52   50.44   51.27   52.79



## SERIE ALLELIQUE ----

### fonctions ----
# Effets série allélique : 

# Pour 1 mk avec une proportion k de variance tot sigma_G expliquée

mc_allele_effect_series <- function(k, nc, sigma_G, denom, s){
  
  f <- 1/2
  u <- sample(c(1, -1), size = 1)
  
  # l'effet 'initial'
  beta1_i <- u * f * sqrt((k* sigma_G)/ denom)
  
  # la liste des facteurs s selon le nombre de croisement :
  list_fact_s<- vector(mode = "list", length = nc)
  list_fact_s[[1]]<-1
  for (r in 2:nc){
    facteur_s <-s^(r-1)
    list_fact_s[[r]] <- facteur_s
  }
  
  # mélanger les effets pour qu'ils soient alloués 
  # aléatoirement à 1 croisement
  facteur_s_bis <- sample(list_fact_s)
  
  beta_effets<- vector(mode = "list", length = nc)
  # les effets de la série allélique (dépendent de l'effet initial)
  for (i in 1:nc){
    beta_effets[[i]]<- beta1_i * facteur_s_bis[[i]]
  }
  
  return(beta_effets)
}

betas <- mc_allele_effect_series(k = k_i[1], nc = 2, 
                                 sigma_G = 50, 
                                 denom = denom, 
                                 s = 0.5^2)


# Determine les phenotypes en ayant les effets betas
mc_calcul_pheno <- function(betas, nc, tab_param_cross){
  
  phenos <- vector(mode='list', length = nc)
  
  for (i in 1:nc){
    
    pheno_i <- tab_param_cross[[i]]$geno_i_sel * betas[[i]]
    phenos[[i]] <- pheno_i
  }
  
  phenob <-phenos[[1]]
  for (i in 2:nc){
    phenob <- rbind(phenob, phenos[[i]])
  }
  
  return(phenob)
}

### BOUCLE ----

n_rep <-100 # nb de rep
nb_mk <- 10
n_cross <-8

# Définir les proportions décroissantes k_i
k_i <- 1 / (1:nb_mk)^2
k_i <- k_i / sum(k_i) # Normaliser pour que la somme soit 1

sigma_G <- 50
s <- 0.9^2

tab_somme <-data.frame(somme = rep(NA, n_rep))

for (o in 1:n_rep){
  print(o)
  
  var_s <- rep(NA, n_rep)
  
  for (k in 1:nb_mk){
    
    
    pos_i <- sample(1:ncol(geno), size = 1)
    # pos_i <- k
    
    tab_param_cross <- mc_param_cross(geno, nc = n_cross,pos_i)
    
    denom <- mc_calcul_denom_s(nc = n_cross, 
                            tab_param_cross = tab_param_cross, 
                            s = s)
    
    # determine the QTL effect
    betas <- mc_allele_effect_series(k = k_i[k], nc = n_cross, 
                                     sigma_G = sigma_G, 
                                     denom = denom, 
                                     s = s )
    # calcule des phenotypes
    pheno_serie <- mc_calcul_pheno(betas = betas, 
                                        nc = n_cross, 
                                        tab_param_cross = tab_param_cross)
    
    var_s[k] <- var(pheno_serie)
    var_res<-data.frame(var_s)
    
  }
  
  tab_somme[o,1] <- sum(var_res$var_s, na.rm = TRUE)
  
  
}


## 2 CROSS : 

# pour 100 rep avec 10 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.15   49.10   51.35   51.13   52.99   55.87 

# pour 100 rep avec 100 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 45.65   49.50   51.30   50.84   52.59   55.01 

# pour 100 rep avec 500 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 44.72   48.91   51.46   51.10   53.04   55.57 

## 5 CROSS : 

# pour 100 rep avec 10 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.09   50.19   51.19   50.91   51.93   53.34

# pour 100 rep avec 100 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 46.69   50.13   51.28   51.05   52.00   53.01 

# pour 100 rep avec 500 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.45   50.35   51.29   51.07   51.95   53.53 


## 8 CROSS : 

# pour 100 rep avec 10 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 48.03   49.91   50.53   50.48   51.29   52.67 

# pour 100 rep avec 100 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.75   50.02   50.72   50.61   51.34   52.80 


# pour 100 rep avec 500 mk : 
boxplot(tab_somme$somme)
abline(h=50, col = 'cyan', lty = 2)
summary(tab_somme$somme)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 48.24   50.11   50.71   50.64   51.27   52.70 



