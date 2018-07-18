#' This function computes the least-squares Procrustes distance between
#' each pair of matrices (configurations of landmarks) from the input set
#'
#'
#' @param X The input set of nx3 matrices (objects)
#'
#'
#' @return The LS Procrustes distance matrix between pairs of objects
#'
#' @author Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
#'
#' @usage
#' cmdistance_RPS(X)
#'
#' @examples
#' source = array(matrix(nrow = 8,ncol = 3),c(8,3,3),dimnames = NULL)
#' source[,,1] <- matrix(c(3,0,0,3,0,1,3,1,1,3,1,0,0,0,0,0,0,1,0,1,1,0,1,0)
#'                    ,nrow = 8,ncol = 3,byrow = TRUE)
#' source[,,2] <- matrix(c(3, 0 ,0,3, 0, 0.5,3, 1 ,0.75,3 ,1 ,0,0 ,0 ,0,0, 0 ,1,0, 1, 1,0, 1, 0.25)
#'                      ,nrow = 8,ncol = 3,byrow = TRUE)
#' source[,,3] <- matrix(c(5, 2 ,1,3, 0, 1.5,3.4, 1 ,1.75,3 ,1 ,0,0 ,0 ,0,0, 2 ,1,0, 3, 1,0, 1, 0.75)
#'                      ,nrow = 8,ncol = 3,byrow = TRUE)
#' result <- RPS::robgit_RPS(source)
#' RPS::cmdistance_RPS(result)
#'
#' @export
cmdistance_RPS <- function(X) {
    numEjemplares <- ncol(X[1, , ])
    nLandmark <- nrow(X[, , 1])


    D <- matrix(nrow = numEjemplares, ncol = numEjemplares, 0)


    for (i in 1:numEjemplares) {
        j <- (i + 1)
        while (j <= numEjemplares) {
            distanciaTotal <- 0
            for (land in 1:nLandmark) {
                # print(norm(t(as.matrix(X[land,,i]))-t(as.matrix(X[land,,j])))^2)
                distanciaLocal = norm(t(as.matrix(X[land, , i])) - t(as.matrix(X[land, , j])))^2
                distanciaTotal = distanciaTotal + distanciaLocal
            }
            D[i, j] <- distanciaTotal^0.5
            D[j, i] <- distanciaTotal^0.5
            j <- j + 1
        }
    }
    return(D)
}
