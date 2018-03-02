################
# IncMatQTLMQE #
################

# function to produce different type of QTL incidence matrices

IncMatQTLMQE <- function(x, mppData, mppData_bi, Q.eff, par.clu,
                               cross.mat, par.mat, order.MAF){

  par.clu <- par.clu[, mppData$parents] # order column of par.clu

  if(Q.eff == "biall") {


    IncMat_QTL(x = x, mppData = mppData_bi, Q.eff = Q.eff, par.clu = par.clu,
                                 cross.mat = cross.mat, par.mat = par.mat)

  } else {

    IncMat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, par.clu = par.clu,
               cross.mat = cross.mat, par.mat = par.mat, order.MAF = order.MAF)

  }

}
