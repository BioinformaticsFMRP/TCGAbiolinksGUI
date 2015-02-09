#' @title Return list with all divisors of the number 
#' @description \code{prime} The function return a array list with all divisors of the number
#' @param x A number
#' @return vector whit \code{x} list of divisors
#' @examples
#' divisor(72)
#' divisor(50)
#' \dontrun{
#'    divisor('a')
#' }
#' @seealso \url{http://en.wikipedia.org/wiki/Divisor}

divisor <- function(number){
  list <- vector()
  for(i in 2:number){
    if((number %% i) == 0){
      list <- c(list, i)
    }
  }
  a <- 1
  return(list)
}