second.round.cbs.matching <- function(diff.data,alpha)
  {
    n <- length(diff.data)
    new.cna <- CNA(diff.data,rep(1,n),1:n)
    new.segment <- segment(new.cna,alpha=alpha)$output
    if(nrow(new.segment)==1) output <- NULL
    else output <- cumsum(new.segment[,5])[-nrow(new.segment)]
    output
  }
