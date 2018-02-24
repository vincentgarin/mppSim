#############
# sign_star #
#############

#' Significance star
#'
#' Function that convert significance of results into stars
#'
#' @param x significance value between 0 and 1.
#'
#' @return Star value
#'
#' @export

sign_star <- function(x){

  if(is.na(x)){

    sign <- NA

  } else {

    if(x>=0.05){
      sign <- ""
    }else if((0.05>x) & (0.01<=x)){
      sign <- "*"
    }else if((0.01>x) & (0.001<=x)){
      sign <- "**"
    }else if(0.001>x){
      sign <- "***"
    }

    return(sign)

  }

}
