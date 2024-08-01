#' LD_computation
#'
#' @description Calculate linkage disequilibrium (LD) between a set of markers
#'
#' @param geno \code{Numeric} genotype marker matrix coded in 0, 1, 2 representing
#' the number of allele copies of one allele (e.g. second parent, minor allele).
#' Genotypes in row and markers in column. The column (marker) name should be
#' identical to the first column of the map (marker indicators).
#' Can be obtrained with \code{\link{sim_mpp_cross}}
#'
#' @param map \code{data.frame} representing the marker map information with
#' the following columns: 1) \code{character} marker identifier;
#' 2) \code{numeric} chromosomes (1, 2, 3, ...);
#' 3) \code{numeric} genetic distance in cM;
#' 4) \code{numeric} Optional physical distance in bp.
#'
#' @param pop_ind Optional population indicator \code{Character vector}
#' indicating to which (sup)population each genotype belongs. Default = NULL
#'
#' @param max_g_dist \code{Numeric value} indicating the maximum genetic distance
#' considered between two markers to calculate the LD. Default = 30
#'
#' @param max_g_dist \code{Numeric value} indicating the selected chromosome.
#' Default = 1
#'
#' @param max_g_dist \code{Numeric value} between 0 and 1 indicating the
#' proportion of marker pairs selected to calculate the LD. Default = 1
#'
#' @return \code{data.frame} with marker i and j position, genetic distance
#' (physical distance) beteen markers and LD
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
#'
#' @export


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
