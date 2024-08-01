library(simcross)
library(mppSim)
library(mppR)
library(LDcorSV)
library(ggplot2)

# sub-functions
source("~/WD/R/packages/mppSim/R/het_det.R")

scale_center_X <- function(X){

  p <- apply(X, 2, FUN = function(x) sum(x)/(2*length(x)))
  P <- rep(1, nrow(X)) %*% (2*t(p))
  Z <- X - P
  den <- sqrt(2 * p * (1 - p))
  Z <- t(t(Z) / den)
  # equivalent to Z %*% diag(1 / den)

  return(Z)

}

# test sim_mpp_cross ----

data("geno_par")
n_par <- 9
rownames(geno_par) <- paste0('P', 1:n_par)

data("EUNAM_map")
rownames(map) <- map[, 1]

cross_scheme_NAM <- function(np){

  n_cross <- np -1

  cr_index <- paste0("%0", nchar(n_cross), "d")
  crosses <- paste0("cr", sprintf(cr_index, 1:n_cross))

  P1 <- rep("P1", n_cross)
  P2 <- paste0("P", 2:np)

  par_per_cross <- cbind(crosses, P1, P2)

  return(par_per_cross)

}

# NAM crossing scheme with 9 parents
cross_scheme <- cross_scheme_NAM(n_par)

# 100 ind per crosses
n_ind_cr_ref <- 100
n_ind_cr <- rep(n_ind_cr_ref, nrow(cross_scheme))

# geno_off <- sim_mpp_cross(geno_par = geno_par,  map = map,
#                           cross_scheme = cross_scheme, n_ind_cr = n_ind_cr)

# open the function

geno_par = geno_par
map = map
cross_scheme = cross_scheme
n_ind_cr = n_ind_cr
geno_score_format = "numeric_IBD"

# possibility to create a geno_score_format argument that specify
# in which format we want to return the genotype scores

# parents_IBS: convert from parents
# parents_IBD: AA BB H (ready to be filled for calc.genoprob)
# numeric_IBD: 1, 2, 3 (parent 1, parent 2, heterozygous within cross meaning)


sim_mpp_cross <- function(geno_par, map, cross_scheme, n_ind_cr,
                          geno_score_format = "parents_IBS"){

  # iteration over the number of crosses

  pop_res <- c()

  n_cross <- nrow(cross_scheme)
  crosses_id <- cross_scheme[, 1]

  map_pos <- split(x = map[, 3], f = factor(map[, 2]))

  n_chr <- length(map_pos)
  chr_lgh <- unlist(lapply(map_pos, max))

  for(i in 1:n_cross){

    # select the two parents scores
    cr_par <- cross_scheme[i, 2:3]
    geno_p1p2 <- geno_par[cr_par, ]

    # determine the heterozygous

    het_score <- apply(X = geno_p1p2, MARGIN = 2, FUN = het_det)
    geno_pipj <- rbind(geno_p1p2[1, ], het_score, geno_p1p2[2, ])

    # simulate the cross genotypes

    f2 <- sim_ril_pedigree(ngen = 1, selfing = TRUE, parents = 1:2)
    ext_ind <- data.frame(id = 5:(n_ind_cr[i]+4), mom=3, dad=3, sex=0,gen=2)
    f2_ext <- rbind(f2, ext_ind)

    xodat_ext <- sim_from_pedigree(f2_ext, L=chr_lgh, xchr = rep(FALSE, n_chr))
    geno_sim <- convert2geno(xodat_ext, map = map_pos)

    # fill with the actual parent genotypes scores

    geno_res_i <- c()

    for(j in 1:n_chr){

      geno_sim_chr_j <- geno_sim[[j]][5:(n_ind_cr[i]+4), ]

      if(geno_score_format == "parents_IBS"){

        for(k in 1:dim(geno_sim_chr_j)[2]){

          geno_sim_chr_j[, k] <- geno_pipj[as.numeric(geno_sim_chr_j[, k]), k]

        }

      } else if (geno_score_format == "numeric_IBD"){

        # geno_sim_chr_j <- apply(geno_sim_chr_j, 2, as.character)
        # geno_sim_chr_j[geno_sim_chr_j == "1"] <- "0"
        # geno_sim_chr_j[geno_sim_chr_j == "3"] <- "1"
        # geno_sim_chr_j <- apply(geno_sim_chr_j, 2, as.numeric)
        geno_sim_chr_j <- geno_sim_chr_j - 1

      }

      geno_res_i <- cbind(geno_res_i, geno_sim_chr_j)

    }

    # add labels

    geno_n_id <- paste0("%0", nchar(n_ind_cr[i]),"d")
    geno_id <- paste0("g", sprintf(geno_n_id, 1:n_ind_cr[i]))
    nm_geno_i <- paste0(crosses_id[i], '_', geno_id)

    rownames(geno_res_i) <- nm_geno_i

    pop_res <- rbind(pop_res, geno_res_i)

  }

  colnames(pop_res) <- map[, 1]
  return(pop_res)

}

geno <- sim_mpp_cross(geno_par = geno_par, map = map,
                      cross_scheme = cross_scheme,
                      n_ind_cr = n_ind_cr,
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
boxplot(res)

abline(h = 12.5, col = "green")
abline(h = 25, col = "blue")

# the divergences are due to the marker that do not segregate according to HWE.

# LD pattern check ----

# get the upper triangle matrix for all pairs of markers
# with genetic and physical distances

# @param geno \code{Numeric} genotype marker matrix coded in 0, 1, 2
# representing the number of allele copies of one allele
# (e.g. second parent, minor allele). Genotypes in row and markers in column.
# The column (marker) name should be identical to the first column of the map
# (marker indicators).

# @param map \code{data.frame} representing the marker map information with
# the following columns: 1) \code{character} marker identifier;
# 2) \code{numeric} chromosomes (1, 2, 3, ...);
# 3) \code{numeric} genetic distance in cM;
# 4) \code{numeric} Optional physical distance in bp.

# @param pop_ind Optional population indicator \code{Character vector}
# indicating to which (sup)population each genotype belongs. Default = NULL

# @param max_g_dist \code{Numeric value} indicating the maximum genetic distance
# considered between two markers to calculate the LD. Default = 30

# @param max_g_dist \code{Numeric value} indicating the selected chromosome.
# Default = 1

# @param max_g_dist \code{Numeric value} between 0 and 1 indicating the
# proportion of marker pairs selected to calculate the LD. Default = 1

# process the data

# remove the markers that are monomorphic
MAF <- apply(geno, 2, FUN = function(x) sum(x)/(2*length(x)))
rem_mk <- MAF == 0 | MAF == 1
map <- map[!rem_mk, ]
geno <- geno[, !rem_mk]

# select (two) crosses
cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)
geno <- geno[cr_ind %in% unique(cr_ind)[1:2], ]
cr_ind <- cr_ind[cr_ind %in% unique(cr_ind)[1:2]]

# pop_ind = cr_ind
# max_g_dist = 30
# chr_sel = 1
# prop_mk_comb = 1


LD_computation <- function(geno, map, pop_ind = NULL, max_g_dist = 30,
                     chr_sel = 1, prop_mk_comb = 1){

  # check data format ----

  if(!identical(colnames(geno), map[, 1])){

    stop("The marker identifiers in geno (column names) and in map (first column) are not identical.")

  }

  if(ncol(map) == 4){ph_dist <- TRUE} else {ph_dist <- FALSE}

  # select chromosome ----
  map <- map[map[, 2] == chr_sel, ]
  geno <- geno[, map[, 1]]

  # get the genetic distances ----

  # marker pair coefficients
  n_mk <- nrow(map)
  m <- matrix(NA, n_mk, n_mk)
  ok <- lower.tri(m, diag = FALSE)
  mk_coeff <- data.frame(row = col(m)[ok], col = row(m)[ok])
  rm(m, ok)

  # get the distances between all pairs of markers
  g_dist <- vector(mode = "list", length = n_mk)
  if(ph_dist){ p_dist <- vector(mode = "list", length = n_mk) }

  for(i in 1:n_mk){
    g_dist[[i]] <- abs(map[i, 3] - map[(i+1):n_mk, 3])
    if(ph_dist){p_dist[[i]] <- abs(map[i, 4] - map[(i+1):n_mk, 4])}
  }

  # convert to vector and remove the last
  g_dist <- unlist(g_dist[-n_mk])
  if(ph_dist){  p_dist <- unlist(p_dist[-n_mk]) }

  if(!ph_dist){
  d_dist <- data.frame(mk_coeff, g_dist)
  } else {
    d_dist <- data.frame(mk_coeff, g_dist, p_dist)
  }

  # keep only the combination of markers with less than 20 cM distances
  d_dist <- d_dist[d_dist[, 3] < max_g_dist, ]

  # restrict by calculating only for every 10th pair of mk
  stp <- round(1 / prop_mk_comb)
  mkp_sel <- seq(1, nrow(d_dist), stp)
  if(length(mkp_sel) == 1){
    stop("The proportion of selected combination 'prop_mk_comb' is too low and lead to the selection of a single marker combination.")
  }
  d_dist <- d_dist[mkp_sel, ]

  # Calculate an optional population structure matrix ----
  if(!is.null(pop_ind)){
    S_mat <- model.matrix(~ -1 + pop_ind)
    # Need to colect (P - 1) population information
    S_mat <- S_mat[, -ncol(S_mat), drop = FALSE]
  }

  # calculate the LD for each selected pairs of markers
  d_dist$LD <- NA

  # takes some time to calculate
  for(i in 1:nrow(d_dist)){

    mk1_pos <- d_dist$row[i]
    mk2_pos <- d_dist$col[i]

    if(!is.null(pop_ind)){

      d_dist$LD[i] <- Measure.R2S(geno[, c(mk1_pos, mk2_pos)],
                                  struc = S_mat, na.presence = FALSE)

    } else {

      d_dist$LD[i] <- Measure.R2(geno[, c(mk1_pos, mk2_pos)],
                                  na.presence = FALSE)

    }

  }

  if(ph_dist){
    colnames(d_dist) <- c("m_i", "m_j", "gen_dist", "ph_dist", "LD")
  } else {
    colnames(d_dist) <- c("m_i", "m_j", "gen_dist", "LD")
  }

  LD_res <- d_dist

  return(LD_res)

}


# test chr 1
LD_chr1 <- LD_computation(geno = geno, map = map, pop_ind = cr_ind,
                    max_g_dist = 40, chr_sel = 1, prop_mk_comb = 1)

# LD_plot

# @param LD_data \code{data.frame} obtained from function \link{\code{LD_computation}}

# @param smooth_curve \code{Logical} value specifying if a smooth curve should be added. Default = TRUE

LD_data = LD_chr1

LD_plot <- function(LD_data){

  # plot with the smooth line
    mod <- loess(LD ~ gen_dist, data=LD_data)
    x_new <- seq(0.01, max(LD_data$gen_dist), 0.01)
    LD_pred <- predict(mod, newdata = data.frame(gen_dist = x_new))
    d_curve <- data.frame(gen_dist = x_new, LD = LD_pred)

    p <- ggplot(data = LD_data, aes(x = gen_dist, y = LD)) +
      geom_point() +
      geom_line(data = d_curve, aes(x = gen_dist, y = LD), col = "red")

    return(p)

  }

LD_plot(LD_chr1)

# LD_decay

threshold = 0.1

LD_decay <- function(LD_data, threshold = 0.1){

  mod <- loess(LD ~ gen_dist, data=LD_data)
  x_new <- seq(0.01, max(LD_data$gen_dist), 0.01)
  LD_pred <- predict(mod, newdata = data.frame(gen_dist = x_new))
  d_curve <- data.frame(gen_dist = x_new, LD = LD_pred)

  LD_dec <- c()

  for(i in 1:length(threshold)){
    cM_i <- d_curve$gen_dist[which(d_curve$LD < threshold[i])[1]]
    LD_dec <- rbind(data.frame(thre_LD = threshold[i], cM = cM_i))
  }

  return(LD_dec)

}
