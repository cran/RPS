#'  Given a n x n distance matrix D (not necessarily Euclidean) and a initial
#'  set X0 (that is, a n x k matrix) of n seeds in k dim, this function finds a set of
#'  n points in k dimensions X (that is, a k x n matrix) using a resistant criterion
#'  such that the n x n matrix Dk of euclidean distances among these new points X
#'  is as close as possible to D.
#'
#' @param D distance matrix n x n to be approximated
#' @param k dimension of output results
#'
#' @return X  A set of n points in k dimensions
#'
#' @author Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
#'
#' @usage
#' resunivMDS_RPS(D,k)
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
#' RPS::resunivMDS_RPS(distance,2)
#'
#' @export
resunivMDS_RPS <- function(D, k = 2) {
    iteraciones <- 10
    tol <- 1e-09
    nl <- nrow(D)  #numer of specimens
    X <- t(randomMatrix(nl, k))
    Dk <- distAllPairsL2(X)
    c <- sum((D - Dk))
    cant <- 0

    for (iter in 1:iteraciones) {
        print(iter)
        for (ii in 1:nl) {
            for (it in 1:floor(sqrt(iter))) {
                Z <- computeIntersections(X, ii, D)
                b <- spatialmed_landmark(t(Z))
                a <- t(b)

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
        c <- sum(sum((D - Dk)))
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

spatialmed_landmark <- function(X) {
  tol = 1e-09

  n <- nrow(X)
  p <- ncol(X)

  m <- matrix(nrow = 1, ncol = p, 0)
  A <- matrix(nrow = n, ncol = p, 0)

  w <- matrix(nrow = 1, ncol = n, 1)
  s <- matrix(nrow = 1, ncol = n, 0)
  aux <- matrix(nrow = 1, ncol = p, 0)
  auxant <- matrix(nrow = 1, ncol = p, 0)

  tdemu <- matrix(nrow = 1, ncol = p, 0)
  rdemu <- matrix(nrow = 1, ncol = p, 0)
  gamagama <- 0
  sensor <- 0

  A <- X


  aux <- apply(A, 2, mean)  #aux <- w%*%A/sum(w)

  h <- 1
  # print(A)

  while ((median(abs(aux - auxant)) > tol) & (h <= 1000)) {
    # print('imprimo s:') print(aux)
    # print('******************************************************') print('AUX: ') print(aux)
    for (k in 1:n) {

      # print('A[k,]') print(t(as.matrix(A[k,]))) print(norm(aux - t(as.matrix(A[k,])) ) )
      if (norm(aux - t(as.matrix(A[k, ]))) == 0) {
        s[1, k] <- tol
        sensor <- 1
      } else {
        s[1, k] <- w[1, k]/norm(aux - t(as.matrix(A[k, ])), "F")
      }

    }
    auxant <- aux

    tdemu <- s %*% A/sum(s)
    rdemu <- s %*% A
    gamagama <- min(1, sensor, norm(rdemu, "F"))
    # print('imprimo s luego del ciclo:') print(s)
    aux <- (1 - gamagama) * tdemu + as.double(gamagama * aux)
    # print('Aux luego de la multiplicacion: ') print(aux)
    h <- (h + 1)

  }
  m <- aux

  return(m)
}

