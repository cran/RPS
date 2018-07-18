#' This function computes the a resistant distance between
#' each pair of matrices from the input set
#'
#' @param X The input set of nx3 matrices (objects)
#'
#' @return This function computes the sum of non-squared euclidean distances across
#'         landmarks for each pair of matrices from the input set
#'
#' @author Guillermo A. Pacheco, Viviana Ferraggine, Sebastian Torcida
#' @usage
#' resdistance_RPS(X)
#'
#'
#' @examples
#' source = array(matrix(nrow = 8,ncol = 3),c(8,3,3),dimnames = NULL)
#' source[,,1] <- matrix(c(3,0,0,3,0,1,3,1,1,3,1,0,0,0,0,0,0,1,0,1,1,0,1,0)
#'                    ,nrow = 8,ncol = 3,byrow = TRUE)
#' source[,,2] <- matrix(c(3, 0 ,0,3, 0, 0.5,3, 1 ,0.75,3 ,1 ,0,0 ,0 ,0,0, 0 ,1,0, 1, 1,0, 1, 0.25)
#'                      ,nrow = 8,ncol = 3,byrow = TRUE)
#' source[,,3] <- matrix(c(5, 2 ,1,3, 0, 1.5,3.4, 1 ,1.75,3 ,1 ,0,0 ,0 ,0,0, 2 ,1,0, 3, 1,0, 1, 0.75)
#'                      ,nrow = 8,ncol = 3,byrow = TRUE)
#' result <- RPS::robgit_RPS(source, consenso = FALSE)
#' RPS::resdistance_RPS(result)
#'
#'
#' @export
resdistance_RPS <- function(X) {
    numEjemplares <- ncol(X[1, , ])
    nLandmark <- nrow(X[, , 1])

    D <- matrix(nrow = numEjemplares, ncol = numEjemplares, 0)

    for (i in 1:numEjemplares) {
        j <- (i + 1)
        while (j <= numEjemplares) {
            distanciaTotal <- 0
            for (land in 1:nLandmark) {
                # print(norm(t(as.matrix(X[land,,i]))-t(as.matrix(X[land,,j]))))
                distanciaLocal = norm(t(as.matrix(X[land, , i])) - t(as.matrix(X[land, , j])))
                distanciaTotal = distanciaTotal + distanciaLocal
            }
            D[i, j] <- distanciaTotal
            D[j, i] <- distanciaTotal
            j <- j + 1
        }
    }
    return(D)
}
