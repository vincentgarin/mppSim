##############
# plot_TP_FP #
##############

#' QTL profile plot with true and false positive QTLs
#'
#' @param Qprof Object of class \code{QTLprof} returned by the function
#' \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}}, or three columns
#' \code{numeric matrix} with chromosome, marker position and -log10(p-values).
#'
#' @param QTL_TRUE Optional argument. List of QTL true QTL positions
#'
#' @param QTL_TP Optional argument. List of false positive QTLs.
#'
#' @param QTL_FP Optional argument. List of QTL true positive QTLs.
#'
#' @param type \code{Character} expression indicating the type of plot should be
#' drawn: "l" for lines , "h" for vertical bar. Default = "l".
#'
#' @param main Title of the graph. Default = "QTL profile".
#'
#' @param threshold \code{Numeric} QTL significance threshold value draw on
#' the plot. Default = 3.
#'
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 18.
#'
#' @author Vincent Garin
#'
#'
#' @examples
#'
#' Not now
#'
#' @export
#'

plot_TP_FP <- function(Qprof, QTL_TRUE, QTL_TP = NULL, QTL_FP = NULL,
                       type = "l", main = "QTL profile", threshold = 3,
                       text.size = 18)
{

  if(inherits(Qprof, "QTLprof")){

    chr <- Qprof$chr
    pos.cM <- Qprof$pos.cM
    log10pval <- Qprof$log10pval
    Qprof <- data.frame(chr, pos.cM, log10pval)

  } else {

    if(!((is.matrix(Qprof)) & (dim(Qprof)[2] == 3) & (is.numeric(Qprof)))){

      stop(paste("The Qprof argument must be a three column numeric matrix with",
                 "chromosome, marker position and -log10(p-value)."))

    }

    Qprof <- data.frame(Qprof)
    colnames(Qprof) <- c("chr", "pos.cM", "log10pval")

  }

  # QTL positions

  pos_Q_TRUE <- QTL_TRUE[, 2:3, drop = FALSE]

  if(!is.null(QTL_TP)){ pos_Q_TP <- QTL_TP[, c(2, 4), drop = FALSE] }

  if(!is.null(QTL_FP)){ pos_Q_FP <- QTL_FP[, c(2, 4), drop = FALSE] }

  if (type == "l") {

    plot <- ggplot(Qprof, aes(x = pos.cM, y = log10pval, group = chr)) +
      geom_line() + facet_wrap(nrow = 1, ~ chr, scales = "free_x") +

      geom_vline(aes(xintercept = pos.cM), pos_Q_TRUE, colour = "black") +

      geom_hline(yintercept = threshold, colour = "red") + theme_bw() +
      xlab("position [cM]") + ylab("-log10(p.val)") +
      ggtitle(main) +

      theme(axis.title.x = element_text(size=text.size),
            axis.title.y = element_text(size=text.size),
            axis.text.x  = element_text(size=text.size - 8),
            axis.text.y = element_text(size = text.size),
            plot.title = element_text(size=(text.size + 4)),
            strip.text.x =  element_text(size=text.size))

    if(!is.null(QTL_TP)){

      plot <- plot + geom_vline(aes(xintercept = pos.cM), pos_Q_TP,
                         linetype = "longdash", colour = "blue")
    }

    if(!is.null(QTL_FP)){

      plot <- plot + geom_vline(aes(xintercept = pos.cM), pos_Q_FP,
                                linetype = "longdash", colour = "red")
    }

  } else if (type == "h") {

    plot <- ggplot(Qprof, aes(x = pos.cM, xend = pos.cM, y = 0, yend = log10pval,
                      group = chr)) + geom_segment() +
      facet_wrap(nrow = 1, ~ chr, scales = "free_x") +

      geom_vline(aes(xintercept = pos.cM), pos_Q_TRUE, colour = "black") +

      geom_hline(yintercept = threshold, colour = "red") +
      theme_bw() + xlab("position [cM]") + ylab("-log10(p.val)") +
      ggtitle(main) +
      theme(axis.title.x = element_text(size=text.size),
            axis.title.y = element_text(size=text.size),
            axis.text.x  = element_text(size=text.size - 8),
            axis.text.y = element_text(size = text.size),
            plot.title = element_text(size=(text.size + 4)),
            strip.text.x =  element_text(size=text.size))

    if(!is.null(QTL_TP)){

      plot <- plot + geom_vline(aes(xintercept = pos.cM), pos_Q_TP,
                                linetype = "longdash", colour = "blue")
    }

    if(!is.null(QTL_FP)){

      plot <- plot + geom_vline(aes(xintercept = pos.cM), pos_Q_FP,
                                linetype = "longdash", colour = "red")
    }

  }

  return(plot)

}
