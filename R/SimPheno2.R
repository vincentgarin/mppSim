############
# SimPheno2 #
############

#' Simulation of phenotypic values
#'
#' Function simulating phenotypic values based on polygenic, QTL and
#' (environmental) error component. The simulation admit 8 types of QTLs: Q1
#' has a different allele in each cross; Q2 segregates in half of the crosses
#' with a different allelic value; Q3 has a different allele for each parent;
#' Q4 has a unique allele hold by a single parent; Q5 has a different allele for
#' each ancestral group; Q6 has a unique allele hold by the most frequent
#' ancestral group; Q7 and Q8 are bi-allelic QTLs attached to a SNP marker.
#'
#' @param QTL list of QTL position obtained with function \code{\link{QposSelect}}
#'
#' @param k Factor representing the ration between QTL (VQ) and polygenic (Vg)
#' variance (VQ = k * Vg). Default = 1.
#'
#' @param her heritability.
#'
#' @param mppData \code{mppData} object
#'
#' @param cross_frac Two values vector indicating the fraction of cross
#' with non-zero value for the cross-specific QTLs Q1 and Q2.
#'
#' @return List with the simulated phenotypic values, each contribution
#' (polygenic, QTL, error), and the proportions of polygenic, QTL, error,
#' QTL covariance, covariance between QTL and polygenic, etc. And also a list
#' of extended QTLs effect with the proportion of non-zero QTL effect
#' contribution (~ QTL allele frequency estimation).
#'
#' @author Vincent Garin
#'
#' @export


SimPheno2 <- function(QTL, k = 1, her, mppData, cross_frac = c(0.5, 0.33)){

  NQTL <- dim(QTL)[1]

  # QTL incidence matrices

  Q.pos <- QTL[, 1]
  Q.pos <- which(mppData$map[, 1] %in% Q.pos)
  Q.eff <- QTL[, 4]

  QTL_inc <- mapply(FUN = inc_mat_QTL, x = Q.pos, Q.eff = Q.eff,
                    MoreArgs = list(mppData = mppData), SIMPLIFY = FALSE)

  # Form the different Beta

  Beta_list <- BetaVal2(QTL = QTL, QTL_inc = QTL_inc, cross_frac = cross_frac)

  # Compute the unscaled realized variance

  vQuS <- mapply(FUN = function(x, y) var(x %*% y), x = QTL_inc, y = Beta_list)

  # obtain the scaling factor

  f <- (QTL[, 6]/sum(QTL[, 6]))/vQuS

  Beta_scaled <- mapply(FUN = function(x, f) x * sqrt(f), x = Beta_list, f = f,
                        SIMPLIFY = FALSE)

  # Recalculate the scaled variance

  vQS <- mapply(FUN = function(x, y) var(x %*% y), x = QTL_inc, y = Beta_scaled)

  # form the QTL contribution per QTL

  y_QTL_i <- mapply(FUN = function(x, y) x %*% y, x = QTL_inc, y = Beta_scaled)

  # check the independence between simulated QTLs contributions

  V_mat <- cov(y_QTL_i)
  diag(V_mat) <- rep(0, NQTL)

  cov_QTL <- sum(c(V_mat))

  y_QTL <- rowSums(y_QTL_i)

  v_QTL <- var(y_QTL)

  v_err <- ((1 - her)/her) * v_QTL

  y_err <- rnorm(n = length(y_QTL), mean = 0, sd = sqrt(v_err))

  Qe_mat <- cbind(y_QTL, y_err)

  cov_Q_e <- 2*cov(Qe_mat)[1, 2]

  ### simulated phenotypic values

  y_sim <- y_QTL + y_err

  # variance components

  v_tot <- var(y_sim)
  v_QTL <- var(y_QTL)
  v_QTLs <- apply(X = y_QTL_i, MARGIN = 2, FUN = var)
  v_err <- var(y_err)

  # proportions

  p_QTL <- v_QTL/v_tot * 100
  p_err <- v_err/v_tot * 100

  v_rest <- cov_QTL + cov_Q_e

  p_rest <- v_rest/v_tot * 100
  p_cov_Q <- cov_QTL/v_tot * 100
  p_cov_Qe <- cov_Q_e/v_tot * 100
  p_Q_6 <- max(v_QTLs)/v_tot * 100
  p_Q_2 <- min(v_QTLs)/v_tot * 100


  res <- list(y_sim = y_sim, y_QTL = y_QTL, y_err = y_err,
              p_QTL = p_QTL, p_err = p_err, p_rest = p_rest,
              p_cov_Q = p_cov_Q, p_cov_Qe = p_cov_Qe, p_Q_6 = p_Q_6,
              p_Q_2 = p_Q_2)

  return(res)

}
