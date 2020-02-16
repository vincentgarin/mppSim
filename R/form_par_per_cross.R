######################
# form_par_per_cross #
######################

#' form parent per cross information in a diallel
#'
#'
#' @param np Parent number in the diallel.
#'
#' @return parent per cross information
#'
#' @author Vincent Garin
#'
#' @export


form_par_per_cross <- function(np = 9){

  # Reference diallel design par per cross

  n_cross <- (np * (np - 1))/2

  crosses <- paste0("cr", sprintf("%02d", 1:n_cross))

  P1 <- paste0("P" , rep(x = 1:(np-1), times = (np-1):1))

  P2_seq <- c()

  for(i in 1:(np - 1)){

    P2_seq_i <- (1 + i):np

    P2_seq <- c(P2_seq, P2_seq_i)

  }

  P2 <- paste0("P", P2_seq)


  par_per_cross <- cbind(crosses, P1, P2)

  return(par_per_cross)

}
