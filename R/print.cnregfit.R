#' Print method for class cnregfit
#' 
#' @param x object of class cnregfit
#' 
#' @export print.cnregfit
#' @export

print.cnregfit = function(x){
  if(!is.null(x$error)) cat("\n");cat(x$error);cat("\n")
  
  object_to_print = 
    data.frame()
  
}

