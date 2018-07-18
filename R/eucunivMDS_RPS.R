#'  Given a n x n distance matrix D (not necessarily Euclidean) and a initial
#'  set X0 of n seeds in k dim (that is, an initial n x k matrix), this function finds
#'  a set of n points in k dimensions X (a final n x k matrix) through a least-squares
#'  criterion such that the n x n matrix Dk of euclidean distances among these new points X
#'  is as close as possible to D.
#'
#' @param D distance matrix n x n to be approximated
#' @param k dimension of output results
#'
#' @return X  A set of n points in k dimensions
#'
#' @author Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
#'
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
#' distance <- RPS::resdistance_RPS(result)
#' RPS::eucunivMDS_RPS(distance,2)
#'
#' @export
eucunivMDS_RPS <- function(D, k = 2) {
    iteraciones <- 10
    tol <- 1e-09
    nl <- nrow(D)
    X <- t(randomMatrix(nl, k))

    Dk <- distAllPairsL2(X)
    c <- sum((D - Dk)^2)
    cant <- 0

    for (iter in 1:iteraciones) {
        for (ii in 1:nl) {
            for (it in 1:floor(sqrt(iter))) {
                print(iter)
                Z <- computeIntersections(X, ii, D)
                a <- t(t(apply(Z, 1, mean)))
                x1 <- cbind(t(t(X[, 1:ii - 1])), a)
                if ((ii + 1) <= nl) {

                  X <- cbind(x1, t(t(X[, (ii + 1):nl])))

                } else {
                  X <- x1
                }

            }
        }


        Dk <- distAllPairsL2(X)
        cant <- c
        c <- sum((D - Dk)^2)
        if (abs(c - cant) < tol) {
            break
        }
    }
    return(t(X))
}

computeIntersections <- function(X, ii, D) {
  k <- nrow(X)
  nl <- ncol(X)
  V <- X - matlab::repmat(t(t(as.matrix(X[, ii]))), 1, nl)
  Daux <- t(t(sqrt((apply(V * V, 2, sum)))))

  Q <- as.matrix(D[ii, ])/Daux  #mirar el operador

  Z <- X - V * matlab::repmat(Q, k, 1)
  Z <- Z[, -ii]
  return(Z)
}

randomMatrix <- function(NRows, NCols) {
  myMat <- matrix(runif(NCols * NRows), ncol = NCols)
  return(myMat)
}

distAllPairsL2 <- function(X) {
  q <- t(X) %*% X
  n <- ncol(q)
  normx <- matlab::repmat(as.matrix(apply(t(X)^2, 1, sum)), 1, n)
  K <- Re(sqrtm(q * (-2) + normx + t(normx)))
  K <- K - (diag(diag(K)))
  return(K)
}

