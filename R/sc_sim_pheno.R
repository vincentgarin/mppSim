#' sc_sim_pheno
#'
#' @description Simulate phenotypes with featured and residual causal variants
#' in a single cross situation. It is possible to add a dilution, with
#' noncausal SNP in the genomic feature.
#'
#' @param X \code{Numeric} genotype marker matrix with genotype x marker
#' information coded as 0, 1, 2 representing the number of copies of the
#' minor allele.
#'
#' @param map genetic marker map: marker id, chromosome, genetic position,
#' physical position.
#'
#' @param n_QTL number of QTLs (causal variants) to be divided into
#' featured and residual variants. Default = 100
#'
#' @param QTL_distribution Fixed to "random" (for the moment)
#'
#' @param min_distance Minimum distance in kb between two selected QTLs.
#' Default = NULL. Not used now
#'
#' @param prop_f Proportion of featured causal variants. Default = 0
#'
#' @param Sp Phenotypic variance. Default = 100
#'
#' @param mu Grand mean. Default = 50
#'
#' @param h2 heritability. Default = 0.5
#'
#' @param h2f Proportion of genetic variance due to the featured (annotated)
#' markers. Default = 0
#'
#' @param scale_mk \code{logical} value specifying if marker scores should be
#' scaled. Default = FALSE
#'
#' @param theta Parameter values of the marker distribution function
#' (\eqn{\theta^{n_{QTL}}}). Default = 0.95
#'
#' @param kj Character string specifying the distribution of proportion
#' by the causal variant in the featured and residual class. One of "constant",
#' "decrease_along_genome", "decrease_random". Default = "decrease_random"
#'
#' @param corr_f \code{numeric} value as correction factor in from of the simulated
#' QTL effects. Default = 1 (no correction correspond to unscaled marker scores).
#'
#' @param dilution_percent \code{numeric} value. Proportion of noncausal markers which are added to
#' the featured markers. Dilution percent must take value between 0 and 100.
#' Default = 0
#'
#' @details
#' The phenotype are simulated following this sequence:
#'
#' 1. Underlying model: \eqn{P_{i} = \mu + G_i + \epsilon_{i}} (unreplicated genotypes) with variance structure
#'
#'\eqn{\sigma_{p}^2 = \sigma_{g}^2 + \sigma_{\epsilon}^2}
#'
#'with \eqn{\sigma_{g}^2 = \sigma_{g_f}^2 + \sigma_{g_r}^2}
#'
#' \eqn{\sigma_{g_f}^2}: genetic variance due to functional SNP
#'
#' \eqn{\sigma_{g_r}^2}: residual genetic variance
#'
#' The latter are linked through the following relationship \eqn{h_f^2 = \frac{\sigma_{g_f}^2}{\sigma_{g}^2} = \frac{\sigma_{g_f}^2}{\sigma_{g_f}^2 + \sigma_{g_r}^2}}
#'
#'  2. Fix \eqn{\sigma_{p}^2 = 100} and \eqn{\mu = 50}
#'
#'  3. Get the genetic variance given assumed variance \eqn{\sigma_{g}^2 = h^2*\sigma_{p}^2}
#'
#'  4. Get the genetic variance under functional annotations \eqn{\sigma_{g_f}^2 = h_f^2*\sigma_{g}^2}
#'
#'  5. Get the residual genetic variance \eqn{\sigma_{g_r}^2 = \sigma_{g}^2 - \sigma_{g_f}^2}
#'
#'  6. Select \eqn{n_{SNP}} that are spread into the \eqn{f} and the \eqn{r} class to form \eqn{X_{f}} and \eqn{X_{r}} given the percentage of markers associated to each category.
#'
#'  7. Calculate for each selected marker the minor allele frequency \eqn{p_j}
#'
#'  8. Determine the marker effect distribution using the following function
#'
#' \eqn{f(\theta, j) = \theta^{j} \quad j=1, ..., n_{SNP}}
#'
#'  9. Determine the proportion of genetic variance accounted for associated to each marker. This is done for each class of effect separately. The values of \eqn{\theta^{j}} are randomly sampled from the general distribution to respect the general distribution of the marker effects.
#'
#'  \eqn{k_j(f) = \frac{\theta^{j}}{\sum_{m=1}^{n_{SNP(f)}}{\theta^{m}}}}
#'
#'  10. For each selected SNP, determine the marker effect to meet the determined \eqn{\sigma_{g_f}^2} and \eqn{\sigma_{g_r}^2} using the formula proposed by Mollandin et al. (2022)
#'
#'  \eqn{\beta_{j(f)} = \pm (\frac{1}{2}) \sqrt{\frac{k_{j(f)}*\sigma_{g_f}^2}{2pq}}}
#'
#'  11. Form the simulated phenotypes
#'
#'  \eqn{\tilde{y} = \mu + X_{s(f)}'\beta_{j(f)} + X_{s(r)}'\beta_{j(r)} + \epsilon_{i}}
#'
#'  where \eqn{X_{s(f)}} (\eqn{X_{s(r)}}) is the centered and scaled marker matrix containing the \eqn{f} (\eqn{r}) markers.
#'
#'  and \eqn{\epsilon_i \sim N(0, \sigma_{\epsilon}^2)} where \eqn{\sigma_{\epsilon}^2 = (1-h^2) * \sigma_{p}^2}
#'
#'  12. Annotation dilution
#'
#'  \eqn{n_{QTL.polluted} = n_{QTL} * \frac{dilution.percent}{100}}
#'
#' Dilution percent :
#' if there is 100 causal markers to be divided into featured and
#' residual variants,
#' 10% of dilution correspond
#' to 10 markers non causal added to the causal featured ones.
#'
#'
#' @return list containing the following object
#' \enumerate{
#' \item {d_y: data.frame with: simulated phenotype (y_sim), featured causal
#' variants contribution (y_f), residual causal variants contribution (y_r),
#' error contribution (e_i)}
#' \item {mk_sel_f: selected featured causal variants
#' (if dilution_percent != 0, Beta = 0 for noncausal variants)}
#' \item {mk_sel_r: selected residual causal variants}
#' \item {X_f: marker matrix of the featured causal variants
#' (if dilution_percent != 0, marker matrix with also noncausal variants
#' featured selected)}
#' \item {X_r: marker matrix of the residual causal variants}
#' \item {Xs_f: standardized marker matrix of the featured causal variants
#' (if dilution_percent != 0, marker matrix with also noncausal variants
#' featured selected)}
#' \item {Xs_r: standardized marker matrix of the residual causal variants}
#' \item {Bf: Additive effect of the featured causal variants}
#' \item {Br: Additive effect of the residual causal variants}
#' }
#' @examples
#'
#' \dontrun{
#' # Generate genotype data
#' data("geno_par")
#' rownames(geno_par) <- paste0('P', 1:nrow(geno_par))
#'
#' data("EUNAM_map")
#' rownames(map) <- map[, 1]
#'
#' # single cross
#' cross_scheme <- cross_scheme_NAM(np = 2)
#' geno_par <- geno_par[1:2, ]
#'
#' geno <- sim_mpp_cross(geno_par = geno_par, map = map,
#'                       cross_scheme = cross_scheme,
#'                       n_ind_cr = rep(100, nrow(cross_scheme)))
#'
#' # without dilution :
#' geno_c1 <- geno$geno_num_IBD
#' y_sim <- sc_sim_pheno(X = geno_c1, map = map)
#'
#' var(y_sim$d_y$y_r)
#'
#' # test simulation consistency
#' res <- matrix(NA, nrow = 100, ncol = 4)
#'
#' for(i in 1:100){
#'
#'   d <- sc_sim_pheno(X = geno_c1, map = map, n_QTL = 1)
#'   res[i, ] <- diag(var(d$d_y))
#'
#'
#' }
#'
#' colnames(res) <- c("Sp", "Sgf", "Sg", "Se")
#' res <- res[, -2]
#' res <- data.frame(res)
#' boxplot(res)
#'
#' abline(h = 100, col = "green")
#' abline(h = 50, col = "blue")
#'
#' # with 10% of dilution (not yet available)
#' # y_sim_dil <- sc_sim_pheno(X = geno, map = map, dilution_percent = 10)
#' }
#'
#' @export

# test values

# data("geno", package = "GPfunctional")
# data("map", package = "GPfunctional")
# X = geno
# map = map
# n_QTL = 100
# QTL_distribution = "random"
# min_distance = 100000
# prop_f = 0.5
# Sp = 100
# mu = 50
# h2 = 0.5
# h2f = 0.5
# theta = 0.95
# kj = "decrease_random"
# dilution_percent=10

sc_sim_pheno <- function(X, map, n_QTL = 100, QTL_distribution = "random",
                          min_distance = NULL, prop_f = 0,
                          Sp = 100, mu = 50, h2 = 0.5, h2f = 0,
                          scale_mk = FALSE,
                          theta = 0.95, kj = "decrease_random",
                          corr_f = 1,
                          dilution_percent = 0){

  # check arguments ----

  if (!is.numeric(dilution_percent)) {
    stop("dilution percent is not a numerical value")
  }
  # if (! ( (dilution_percent<=100) & (dilution_percent>=0) )){
  #   stop("dilution percent must take value between 0 and 100")
  # }

  # get the genetic variance ----
  Sg <- h2 * Sp

  # correct h2f for extreme 0 or 1 proportion
  if(prop_f == 0){h2f <- 0} else if (prop_f == 1){h2f <- 1}

  # partition the genetic variance
  Sgf <- h2f * Sg
  Sgr <- Sg - Sgf


  # sample and partition markers
  map_sel <- QTL_pos_sample(map = map, n_QTL = n_QTL,
                            QTL_distribution = QTL_distribution,
                            min_distance = min_distance)

  # split QTLs into featured and residuals
  map_sel_split <- split_FA_nFA(map_sel = map_sel, prop = prop_f)

  X_f <- X[, map_sel_split$map_FA[, 1], drop = FALSE]
  X_r <- X[, map_sel_split$map_nFA[, 1], drop = FALSE]
  n_mk_f <- ncol(X_f)
  n_mk_r <- ncol(X_r)

  # maps of featured markers and residuals markers
  map_sel_f <- map_sel_split$map_FA
  map_sel_r <- map_sel_split$map_nFA

  # calculate the marker frequency
  p_f <- apply(X_f, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x)))
  p_r <- apply(X_r, MARGIN = 2, FUN = function(x) sum(x)/(2*length(x)))

  # center and scale the marker matrix
  if(scale_mk) Xs <- scale_center_X(X = X) else Xs <- X

  Xs_f <- Xs[, map_sel_f[, 1], drop = FALSE]
  Xs_r <- Xs[, map_sel_r[, 1], drop = FALSE]

  # Determination of the marker dist
  mk_dist <- mk_eff_dist(n = n_QTL, theta = theta)

  if(prop_f > 0){
    sel_kf <- sort(sample(x = 1:n_QTL, size = n_mk_f))
    sel_kr <- -sel_kf
  } else{
    sel_kf <- NULL
    sel_kr <- 1:n_QTL
  }


  if(kj == "constant"){

    k_f <- rep(1/n_mk_f, n_mk_f)
    k_r <- rep(1/n_mk_r, n_mk_r)

  } else if (kj == "decrease_along_genome"){

    k_f <- mk_dist[sel_kf]/sum(mk_dist[sel_kf])
    k_r <- mk_dist[sel_kr]/sum(mk_dist[sel_kr])

  } else if (kj == "decrease_random"){

    k_f <- sample(mk_dist[sel_kf]/sum(mk_dist[sel_kf]))
    k_r <- sample(mk_dist[sel_kr]/sum(mk_dist[sel_kr]))

  }

  # check dist
  # plot(y = mk_dist/sum(mk_dist), x = 1:n_QTL, main = "n=100, theta=0.95")
  # plot(y = k_f, x = 1:length(k_f), main = "k_f distribution")
  # plot(y = k_r, x = 1:length(k_r), main = "k_r distribution")

  # determine the beta_j(f)
  Bjf <- rep(NA, n_mk_f)

  if(length(Bjf) > 0){

    for(f in 1:n_mk_f){
      Bjf[f] <- allele_effect(p = p_f[f], Sg = Sgf, k = k_f[f], f = corr_f)
    }

    map_sel_f$Beta <- Bjf

    # form the simulated phenotype values
    y_f <- (Xs_f %*% matrix(Bjf, ncol = 1))[, 1]

  } else{

    y_f <- rep(0, nrow(Xs_f))

  }

  # determine the beta_j(r)
  Bjr <- rep(NA, n_mk_r)

  if(length(Bjr) > 0){

    for(r in 1:n_mk_r){
      Bjr[r] <- allele_effect(p = p_r[r], Sg = Sgr, k = k_r[r], f = corr_f)
    }

    map_sel_r$Beta <- Bjr

    # form the simulated phenotype values
    y_r <- (Xs_r %*% matrix(Bjr, ncol = 1))[, 1]

  } else{

    y_r <- rep(0, nrow(Xs_r))

  }


  # get the residual
  Se <- (1-h2) * Sp
  e_i <- rnorm(n = nrow(X), mean = 0, sd = sqrt(Se))

  y_sim <- mu + y_f + y_r + e_i

  ## pollution of annotations = dilution
  # choose markers randomly in the genome (not featured or residuals)
  # add them to the group of featured markers

  if (dilution_percent != 0) {

    # Map of markers not selected
    map_dilution <- map[!map[, 1] %in% map_sel[, 1],]

    # Check if map_dilution[, 1] is not empty
    if (length(map_dilution[, 1]) == 0) {
      stop("Aucun marqueur trouvé dans map_dilution[, 1]")
    }

    # Check if markers in map_dilution are in X
    markers_not_in_X <- setdiff(map_dilution[, 1], colnames(X))
    if (length(markers_not_in_X) > 0) {
      stop(paste("Les marqueurs suivants ne sont pas présents dans X:", paste(markers_not_in_X, collapse = ", ")))
    }

    # Matrix of genotypes for markers not selected
    X_no_sel <- X[, map_dilution[, 1]]

    # Number of markers to select
    nb_dilution <- (dilution_percent * n_QTL) / 100

    # Map of selected markers
    map_sel_dilution <- QTL_pos_sample(map = map_dilution, n_QTL = nb_dilution,
                                       QTL_distribution = QTL_distribution,
                                       min_distance = min_distance)

    # New matrix of featured markers
    X_f_dilution <- X[, map_sel_dilution[, 1]]

    X_f <- cbind(X_f, X_f_dilution)

    # Center and scale the marker matrix
    Xs_f_dilution <- Xs[, map_sel_dilution[, 1]]

    Xs_f <- cbind(Xs_f, Xs_f_dilution)

    # Noncausal SNP have a beta = 0
    map_sel_dilution$Beta <- 0
    map_sel_f <- rbind(map_sel_f, map_sel_dilution)
  }


  # check variance components
  var(y_sim)
  var(y_f)
  var(y_r)
  var(e_i)

  return(list(d_y = data.frame(y_sim, y_f, y_r, e_i),
              mk_sel_f = map_sel_f,
              mk_sel_r = map_sel_r,
              X_f = X_f, X_r = X_r, Xs_f = Xs_f,
              Xs_r = Xs_r, Bf = Bjf, Br = Bjr))

}
