#' This function performs the classical least-squares (LS)
#' Procrustes superimposition of the input configurations of landmarks
#'
#' @param X A s-dimensional array (s=2 or s=3) of n x k matrices, representing shapes of k objects through n landmarks in s dimensions
#'
#' @return s-dimensional array of n x k matrices, representing the different specimens after adjusting.
#'
#' @author Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
#'
#' @usage
#' procrustesCM_RPS(X)
#'
#' @examples
#' source = array(matrix(nrow = 8,ncol = 3),c(8,3,3),dimnames = NULL)
#' source[,,1] <- matrix(c(3,0,0,3,0,1,3,1,1,3,1,0,0,0,0,0,0,1,0,1,1,0,1,0)
#'                    ,nrow = 8,ncol = 3,byrow = TRUE)
#' source[,,2] <- matrix(c(3, 0 ,0,3, 0, 0.5,3, 1 ,0.75,3 ,1 ,0,0 ,0 ,0,0, 0 ,1,0, 1, 1,0, 1, 0.25)
#'                      ,nrow = 8,ncol = 3,byrow = TRUE)
#' source[,,3] <- matrix(c(5, 2 ,1,3, 0, 1.5,3.4, 1 ,1.75,3 ,1 ,0,0 ,0 ,0,0, 2 ,1,0, 3, 1,0, 1, 0.75)
#'                      ,nrow = 8,ncol = 3,byrow = TRUE)
#' result <- RPS::procrustesCM_RPS(source)
#' result
#'
#'
#' @export
procrustesCM_RPS <- function(X) {
  result <- geomorph::gpagen(X)
  return(result['coords'][[1]])
}
