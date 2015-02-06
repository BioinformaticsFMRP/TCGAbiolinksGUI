#' @title Verify if a number is prime 
#' @description \code{prime} verify if a number is prime
#' @param x A number
#' @return TRUE if \code{x} is prime, FALSE otherwise
#' @examples
#' prime(4)
#' prime(13)
#' \dontrun{
#'    prime('a')
#' }
#' @seealso \url{http://en.wikipedia.org/wiki/Prime_number}
prime <- function(x){
  if(is.na(x)) {
      stop("x must be a number")
  }
  if(x==2) return (TRUE)
  for(i in 2:sqrt(x)){
    if(x %% i == 0)
      return (FALSE)
  }
  return (TRUE)
}
