###################
# SampleMPPDesign #
###################

#' Sample MPP design
#'
#' Function to sample the realisation of a MPP design from a full diallel made
#' from \code{np_dia} with n_ind individual per crosses. The function sample a
#' realisation of the design specified in argument \code{MPP_des} composed from
#' \code{np_sel} parents including \code{N} individuals.
#'
#' @param np_dia Parent number in the reference diallel. Default = 9.
#'
#' @param n_ind Number of individual per cross in the reference diallel.
#' Default = 200.
#'
#' @param MPP_des Selected MPP design for the MPP sampling. One of: 'NAM' for
#' nested association mapping design; 'Dia' for half-diallel; 'Chess' for
#' chessboard design; 'Fact' for factorial desian; and 'Real' for a realized
#' design.
#'
#' @param np_sel Number of parents in the selected MPP design
#'
#' @param N Total number of sampled individuals. Default = 800.
#'
#' @return List containing the numeric genotype identifiers, the selected parents
#' and the selected crosses.
#'
#' @author Vincent Garin
#'
#' @export

# np_dia <- 9
# n_ind <- 200
# MPP_des <- c("NAM", "Dia", "Chess", "Fact", "Real")
# np_sel <- 9
# N <- 800


SampleMPPDesign <- function(np_dia = 9, n_ind = 200, MPP_des, np_sel, N = 800){

  # Reference diallel design par per cross

  n_cross <- (np_dia * (np_dia - 1))/2

  crosses <- paste0("c", 1:n_cross)

  P1 <- paste0("P" , rep(x = 1:(np_dia-1), times = (np_dia-1):1))

  P2_seq <- c()

  for(i in 1:(np_dia - 1)){

    P2_seq_i <- (1 + i):np_dia

    P2_seq <- c(P2_seq, P2_seq_i)

  }

  P2 <- paste0("P", P2_seq)


  par_per_cross <- cbind(crosses, P1, P2)

  # cross indicator and corresponding numeric indicator

  cross_ind <- rep(crosses, each = n_ind)
  num_ind <- 1:length(cross_ind)

  # Reference identifiers of the crosses

  cr_id <- paste0(par_per_cross[, 2], par_per_cross[, 3])

  if(MPP_des == "NAM"){

    sel_par <- sort(sample(1:np_dia, np_sel))

    # Get the crosses id of the selected design

    cent_par <- sample(sel_par, 1)
    per_par <- sel_par[sel_par != cent_par]

    cr_id_sel <- rep("", (np_sel - 1))

    for(i in 1:(np_sel - 1)){

      par_i <- sort(c(cent_par, per_par[i]))

      cr_id_sel[i] <- paste0("P", paste(par_i, collapse = "P"))

    }

    cr_ind_sel <- par_per_cross[cr_id %in% cr_id_sel, 1]

  } else if (MPP_des == "Dia"){

    sel_par <- sort(sample(1:np_dia, np_sel))
    sel_par_ext <- paste0("P", sel_par)

    cr_id_sel <- outer(X = sel_par_ext, Y = sel_par_ext, FUN = paste0)
    cr_id_sel <- cr_id_sel[upper.tri(cr_id_sel)]

    cr_ind_sel <- par_per_cross[cr_id %in% cr_id_sel, 1]

  } else if (MPP_des == "Chess"){

    sel_par <- sort(sample(1:np_dia, np_sel))
    sel_par_ext <- paste0("P", sel_par)

    cr_id_sel <- outer(X = sel_par_ext, Y = sel_par_ext, FUN = paste0)

    # Generate a checkerboard matrix

    chessMatrix <- solve(Hilbert(np_sel)) < 0
    chessMatrix[lower.tri(chessMatrix)] <- FALSE

    # cr_id_sel <- cr_id_sel[chessMatrix@x] # Problem there
    cr_id_sel <- cr_id_sel[as.matrix(chessMatrix)]

    cr_ind_sel <- par_per_cross[cr_id %in% cr_id_sel, 1]


  } else if (MPP_des == "Fact"){

    sel_par <- sort(sample(1:np_dia, np_sel))

    if(np_sel %% 2 == 0){

      sel_par1 <- paste0("P", sel_par[1:(np_sel/2)])
      sel_par2 <- paste0("P", sel_par[((np_sel/2) + 1):np_sel])

    } else {

      sel_par1 <- paste0("P", sel_par[1:floor(np_sel/2)])
      sel_par2 <- paste0("P", sel_par[ceiling(np_sel/2):np_sel])

    }

    cr_id_sel <- c(outer(X = sel_par1, Y = sel_par2, FUN = paste0))

    cr_ind_sel <- par_per_cross[cr_id %in% cr_id_sel, 1]


  } else if (MPP_des == "Real"){

    sel_par <- 1:9
    cr_ind_sel <- c("c6", "c13", "c19", "c20", "c25", "c28", "c29", "c30",
                    "c31", "c32", "c33")

  }

  # compute the individual per cross

  n_cr <- length(cr_ind_sel)

  if (MPP_des == "Real"){

    ind_cr <- c(24, 48, 73, 99, 109, 104, 115, 70, 40, 39, 30)
    ind_prop <- ind_cr/sum(ind_cr)
    ind_cr <- round(ind_prop * N)
    rest <- sum(ind_cr) - N

    if(rest !=0){

      if(rest < 0) { # Add ind

        id_samp <- sample(1:n_cr, 1)
        ind_cr[id_samp] <- ind_cr[id_samp] + rest

      } else { # Remove ind
        id_samp <- sample(1:n_cr, 1)
        ind_cr[id_samp] <- ind_cr[id_samp] - rest

      }

    }

  } else { # Other MPP designs than the realized one.

    rest <- N %% n_cr

    if(rest == 0){ # the total can be divided roundly into all crosses

      int_div <- N %/% n_cr
      ind_cr <- rep(int_div, n_cr)

    } else {

      int.div <- N %/% n_cr

      ind_cr <- c(rep(int.div, (n_cr-rest)), rep((int.div + 1), (rest)))
      ind_cr <- ind_cr[sample(1:n_cr)]

    }

  }

  # Final indicator sampling given the list of crosses and the numeric indicator

  sel_ind <- c()

  for(i in 1:n_cr){

   num_ind_i <- num_ind[cross_ind %in% cr_ind_sel[i]]
   sel_ind_i <- sample(x = num_ind_i, size = ind_cr[i])

   sel_ind <- c(sel_ind, sel_ind_i)

  }

  sel_ind <- sort(sel_ind)

  return(list(sel_ind = sel_ind, sel_par = sel_par, cr_ind_sel = cr_ind_sel,
              ind_cr = ind_cr))

}
