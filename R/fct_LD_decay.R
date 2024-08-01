#' LD_decay
#'
#' @description Calculate linkage disequilibrium decay
#'
#' @param LD_data \code{data.frame} obtained from function \code{\link{LD_computation}}
#'
#' @param threshold \code{Numerical} value of vector of LD thresholds. Default = 0.1
#'
#' @return \code{data.frame} LD decay limit(s) in genetic distance per LD threshold
#'
#' @examples
#'
#' data("geno_par")
#' rownames(geno_par) <- paste0('P', 1:nrow(geno_par))
#'
#' data("EUNAM_map")
#' rownames(map) <- map[, 1]
#'
#' # NAM crossing scheme with 3 parents
#' cross_scheme <- cross_scheme_NAM(3)
#'
#' geno <- sim_mpp_cross(geno_par = geno_par, map = map,
#'                      cross_scheme = cross_scheme,
#'                      n_ind_cr = c(100, 100),
#'                      geno_score_format = "numeric_IBD")
#'
#' cr_ind <- substr(x = rownames(geno), start = 1, stop = 3)
#'
#' LD_chr1 <- LD_computation(geno = geno, map = map, pop_ind = cr_ind,
#'                          max_g_dist = 40, chr_sel = 1, prop_mk_comb = 1)
#'
#' LD_plot(LD_chr1)
#'
#' LD_decay_01 <- LD_decay(LD_chr1)
#'
#' @export

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
