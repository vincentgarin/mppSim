#' split_FA_nFA
#'
#' @description function to randomly split QTLs into feature and
#' residual class
#'
#' @return list of map (data.frame) per class.
#'
#' @noRd

split_FA_nFA <- function(map_sel, prop){

  if(prop > 0){

    mk_an <- sort(sample(x = 1:nrow(map_sel), size = round(nrow(map_sel)*prop)))
    map_sel_FA <- map_sel[mk_an, ]
    map_sel_nFA <- map_sel[-mk_an, ]

  } else if (prop == 0) {

    mk_Nan <- 1:nrow(map_sel)
    map_sel_FA <- map_sel[-mk_Nan, ]
    map_sel_nFA <- map_sel[mk_Nan, ]

  }

  return(list(map_FA = map_sel_FA, map_nFA = map_sel_nFA))

}
