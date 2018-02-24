################
# QTLDetection #
################

#' QTL detection
#'
#' Function to perform QTL detection using four types of models: cross-specific,
#' parental, ancestral and bi-allelic. The function save the results of the
#' QTL detection processes in a specified folder and return the detected QTLs.
#'
#' @param mppData \code{mppData} object for cross-specific, parental and
#' ancestral model.
#'
#' @param mppData_bi \code{mppData} object for bi-allelic model.
#'
#' @param par_clu Parent clustering object.
#'
#' @param thre_cr Significance thresold for the cross-specific model.
#'
#' @param thre_par Significance threshold for the parental model.
#'
#' @param thre_anc Significance threshold for the ancestral model.
#'
#' @param thre_biall Significance threshold for the bi-allelic model.
#'
#' @param win.cof Minimum distance between two cofactors. Default = 50.
#'
#' @param win.QTL Minimum distance between two selected QTLs. Default = 30.
#'
#' @param folder Folder location where the results are saved.
#'
#' @param Rep_id Identifier for the specific replication.
#'
#' @param parllel If the analysis should be run in parallel
#'
#' @param cluster Cluster object to run the analysis in parallel.
#'
#' @return list with the QTLs detected with the four models.
#'
#' @author Vincent Garin
#'
#' @export

# thre_cr <- 3
# thre_par <- 3.5
# thre_anc <- 4
# thre_biall <- 4.5
# win.cof <- 50
# win.QTL <- 30
# folder <- "results/QTL_detection/MPP_1"
# Rep_id <- "Test"
# parallel <- TRUE
# cluster <- makeCluster(6)

QTLDetection <- function(mppData, mppData_bi, par_clu, thre_cr, thre_par,
                         thre_anc, thre_biall, win.cof = 50, win.QTL = 30,
                         folder, Rep_id, parallel, cluster){

  # Create a folder to store the results

  fold_loc <- file.path(folder, Rep_id)
  # system(paste("mkdir", fold_loc))
  dir.create(fold_loc)

  # adapt the par_clu object

  par_clu_new <- par_clu[, mppData$parents]
  par_clu_new <- par_clu_new[mppData$map[, 1], ]

  par_clu_new <- parent_clusterCheck(par.clu = par_clu_new)[[1]]

  # cross-specific model

  QTL_detect <- tryCatch(mpp_proc(pop.name = "SIM", trait.name = "T",
                                  mppData = mppData, Q.eff = "cr",
                                  thre.cof = thre_cr, thre.QTL = thre_cr,
                                  win.cof = win.cof, win.QTL = win.QTL,
                                  N.cim = 1, output.loc = fold_loc,
                                  parallel = parallel, cluster = cluster,
                                  verbose = FALSE), error = function(e) "error")


  if(is.null(QTL_detect)){ QTL_cr <- "no_QTL"

  } else {

    if(QTL_detect[1] == "error"){ QTL_cr <- "error"; print("error")

    } else { QTL_cr <- QTL_detect$QTL }

  }

  rm(QTL_detect)


  # parental model

  QTL_detect <- tryCatch(mpp_proc(pop.name = "SIM", trait.name = "T",
                                  mppData = mppData, Q.eff = "par",
                                  thre.cof = thre_par, thre.QTL = thre_par,
                                  win.cof = win.cof, win.QTL = win.QTL,
                                  N.cim = 1, output.loc = fold_loc,
                                  parallel = parallel, cluster = cluster,
                                  verbose = FALSE), error = function(e) "error")


  if(is.null(QTL_detect)){ QTL_par <- "no_QTL"

  } else {

    if(QTL_detect[1] == "error"){ QTL_par <- "error"; print("error")

    } else { QTL_par <- QTL_detect$QTL }

  }

  rm(QTL_detect)

  # ancestral model

  QTL_detect <- tryCatch(mpp_proc(pop.name = "SIM", trait.name = "T",
                                  mppData = mppData, Q.eff = "anc",
                                  par.clu = par_clu_new,
                                  thre.cof = thre_anc, thre.QTL = thre_anc,
                                  win.cof = win.cof, win.QTL = win.QTL,
                                  N.cim = 1, output.loc = fold_loc,
                                  parallel = parallel, cluster = cluster,
                                  verbose = FALSE), error = function(e) "error")


  if(is.null(QTL_detect)){ QTL_anc <- "no_QTL"

  } else {

    if(QTL_detect[1] == "error"){ QTL_anc <- "error"; print("error")

    } else { QTL_anc <- QTL_detect$QTL }

  }

  rm(QTL_detect)

  # bi-allelic model

  QTL_detect <- tryCatch(mpp_proc(pop.name = "SIM", trait.name = "T",
                                  mppData = mppData_bi, Q.eff = "biall",
                                  thre.cof = thre_biall, thre.QTL = thre_biall,
                                  win.cof = win.cof, win.QTL = win.QTL,
                                  N.cim = 1, output.loc = fold_loc,
                                  parallel = parallel, cluster = cluster,
                                  verbose = FALSE), error = function(e) "error")

  if(is.null(QTL_detect)){ QTL_biall <- "no_QTL"

  } else {

    if(QTL_detect[1] == "error"){ QTL_biall <- "error"; print("error")

    } else { QTL_biall <- QTL_detect$QTL }

  }

  rm(QTL_detect)

  return(list(QTL_cr = QTL_cr, QTL_par = QTL_par, QTL_anc = QTL_anc,
              QTL_biall = QTL_biall))

}
