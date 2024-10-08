#' QTL_pos_sample
#'
#' @description Sample causal variants (QTL)
#'
#' @return map (data.frame) of the selected variants
#'
#' @noRd

QTL_pos_sample <- function(map, n_QTL, QTL_distribution = "random",
                           min_distance = NULL){

  # copy the map to modify it later
  map_cp <- map

  if(QTL_distribution == "random"){

    mk_sel <- rep(NA, n_QTL)

    for(i in 1:n_QTL){
      # select a position
      p_i <- sample(1:nrow(map_cp), size = 1)
      # add the marker to the selected one
      mk_sel[i] <- map_cp[p_i, 1]

      # discard position on left and right given min distance on the same chr
      # NOT YET
      # map_cp_chr_i <- map_cp[map_cp$chr == map_cp$chr[p_i], ]
      # bp_p_i <- map_cp$phy_pos[p_i]
      # ex_win_i <- abs(bp_p_i - map_cp_chr_i$phy_pos) < min_distance
      # p_ex_i <- map_cp_chr_i$mk.names[ex_win_i]
      # map_cp <- map_cp[!(map_cp$mk.names %in% p_ex_i),]
      p_ex_i <- mk_sel[1:i]
      map_cp <- map_cp[!(map_cp[, 1] %in% mk_sel), ]

      # test if there are still some marker to select
      if(nrow(map_cp) == 0){
        mes_stop <- paste("no more marker to select after", i, "position")
        stop(mes_stop)
      }

    }

    map_sel <- map[map[, 1] %in% mk_sel, ]

  } else if (QTL_distribution == "clustered"){

    # Think about the biological meaning of that scenario
    stop("Not available for the moment.")
  }

  return(map_sel)

}
