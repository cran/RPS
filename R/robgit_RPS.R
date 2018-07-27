#' Resistant Procrustes Superposition Package in R (RPS):
#' a novel package for landmark-based resistant shape analysis
#'
#' RPS provides a set of tools to perform a rather complete descriptive
#' landmark-based resistant shape analysis 3D and 2D, following Torcida
#' et al. 2014 ("An integrated approach for landmark-based resistant shape analysis
#'    in 3D", Evol. Biol. 41(2):351_366). More specifically, these tools enable to
#' obtain: i) a generalized resistant Procrustes superposition (robgit_RPS.R) for
#' a set of configurations of landmarks either in 3D and 2D; ii) a resistant
#' distance (resdistance_RPS.R) to quantify shape differences obtained following
#' the resistant Procrustes superimposition, and iii) a resistant ordination
#' (resunivMDS_RPS.R) of the superimposed configurations based on the universal
#' Multidimensional Scaling from (Agarwal et al. 2010). Corresponding least
#' squares (LS) counterparts of all these tools (procrustesCM_RPS.R,
#' cmdistance_RPS.R and eucunivMDS_RPS.R, respectively) have also been implemented
#' in RPS_R to offer a more complete and self-contained set of shape analysis
#' descriptive tools. This enables the comparison of the LS and resistant
#' superimposition results when  applied to the same dataset.
#' Also included is a rather new method for a
#' resistant analysis of individual shape asymmetry for configurations of
#' landmarks in 2D with bilateral symmetry (matching or object symmetry), following
#' Torcida et al. 2016 ("A resistant method for landmark-based analysis of
#' individual asymmetry in two dimensions",Quant. Biol. 4(4):270_282). The main
#' tools enable to estimate the resistant symmetric shape under matching symmetry
#' (matchingsymm_RPS.R) adn the resistant symmetric shape estimation under object
#' symmetry (objectsymm_RPS.R). In both cases, a plot of the results and the table
#' sof landmarks contributions to asymmetry are also offered.
#'
#' @section Functions:
#' eucunivMDS_RPS, resunivMDS_RPS, cmdistance_RPS, resdistance_RPS, readlandtxtMorphJ_RPS, robgit_RPS,matchingsymm_RPS,objectsymm_RPS, procrustesCM_RPS
#'
#' @docType package
#' @name 'RPS'
#' @param X A s-dimensional array of n x k matrices (k configurations of n landmarks), each representing the shape of an object
#' @param consenso A logical value that determines if the consensus configuration is returned.
#' @return s-dimensional array of n x k matrices, representing the (resistant) superimposed objects
#'
#'
#' @author Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
#'
#' @usage
#' robgit_RPS(X, consenso = FALSE)
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
#' result
#'
#' @export
robgit_RPS <- function(X, consenso = FALSE) {

  adjustament = FALSE
  if (ncol(X[, , 1]) == 2) {
    X <- adjustment3D(X)
    adjustament = TRUE
  }
  tol <- 1e-09
  iterTotal <- 20  #maximum number of iterations
  Z <- X
  f = nrow(Z[, , 1])
  c = ncol(Z[, , 1])
  r = ncol(Z[, 1, ])
  res <- matrix(nrow = f, ncol = r, 0)
  resme <- matrix(nrow = 1, ncol = r, 0)
  restot <- 0
  auxres <- matrix(nrow = f, ncol = r, 0)
  auxresme <- matrix(nrow = 1, ncol = r, 0)
  mejora <- matrix(nrow = f + 1, ncol = r, 0)
  Y <- matrix(nrow = f, ncol = c, 0)
  W <- matrix(nrow = f, ncol = c, 0)
  R <- array(matrix(nrow = f, ncol = c, 0), c(f, c, r + 1), dimnames = NULL)
  X <- array(matrix(nrow = f, ncol = c, 0), c(f, c, r), dimnames = NULL)
  X0 <- X
  E <- matrix(nrow = f, ncol = f, 0)
  A <- matrix(nrow = f, ncol = f, 0)
  Aux <- array(matrix(nrow = f, ncol = c, 0), c(f, c, r), dimnames = NULL)
  H <- array(matrix(nrow = 3, ncol = 3, 0), c(3, 3, r), dimnames = NULL)
  p <- matrix(nrow = 1, ncol = r, 1)
  TR <- matrix(nrow = f, ncol = c, 0)
  X <- Z
  s <- 0
  X <- initialFocus(X, f, r)
  X <- adjustmentScale(X, r)
  Y <- spatialmed_config(X)

  for (k in 1:r) {
    for (i in 1:f) {
      res[i, k] <- norm(t(as.matrix(Y[i, ])) - t(as.matrix(X[i, , k])), "F")
    }
    #resme[1, k] <- median(res[, k])
    resme[1, k] <- sum(res[, k])
  }

  residual <- sum(resme)
  residualant <- 0
  z <- 1
  Aux <- X
  conteo <- matrix(nrow = f + 1, ncol = r, 0)
  conteomed <- matrix(nrow = 1, ncol = r, 0)

  # median landmark residual between consecutive (resistant) consensus is the chosen stopping criterion
  while ((z <= iterTotal) & (abs(residual - residualant) > tol)) {
    for (k in 1:r) {
      Aux[, , k] <- scaleSpecimen(Aux[, , k], Y)
      list_out <- rotation(Aux[, , k], H[, , k], Y)
      Aux[, , k] <- list_out[[1]]
      H[, , k] <- list_out[[2]]
      Aux[, , k] <- translation(Aux[, , k], Y, f)
      for (i in 1:f) {
        auxres[i, k] <- norm(t(as.matrix(Y[i, ])) - t(as.matrix(Aux[i, , k])), "F")
        if (auxres[i, k] <= res[i, k]) {
          conteo[i, k] <- 1
        } else {
          conteo[i, k] <- 0
        }
      }
      auxresme[1, k] <- median(auxres[, k])
      auxresme[1, k]
    }
    for (k in 1:r) {
      if (auxresme[1, k] <= resme[1, k]) {
        X[, , k] <- Aux[, , k]
        conteomed[1, k] <- 1
      } else {
        conteomed[1, k] <- 0
      }
    }
    conteo[f + 1, ] <- matrix(nrow = 1, ncol = r, 0)
    conteo[f + 1, ] <- (100 * apply(conteo, 2, sum))/f
    for (k in 1:r) {
      for (i in 1:f) {
        res[i, k] <- norm(t(as.matrix(Y[i, ])) - t(as.matrix(X[i, , k])), "F")
      }
      #resme[1, k] <- median(res[, k])
      resme[1, k] <- sum(res[, k])
    }

    W <- Y
    Y <- spatialmed_config(X)  # consensus: the spatial median configuration
    residualant <- residual
    residual <- sum(resme)
    z <- (z + 1)

  }

  if (adjustament == TRUE) {
    V <- array(matrix(nrow = f, ncol = 2, 0), c(f, 2, r + 1), dimnames = NULL)
    for (i in 1:r) {
      U <- X[, , i]
      U <- U[, -3]
      V[, , i] <- U
    }
    U <- Y
    U <- U[, -3]
    if (consenso) {
      V[, , r + 1] <- U
      return(V)
    } else {
      return(V[, , -(r + 1)])
    }

  } else {
    R[, , 1:r] <- X
    if (consenso) {
      R[, , r + 1] <- Y
      return(R)
    } else {
      return(R[, , -(r + 1)])
    }

  }
}

adjustment3D <- function(X) {
  f <- nrow(X[, , 1])  # number of landmarks in every object
  r <- ncol(X[, 1, ])  # number of objects in the dataset

  R <- array(matrix(nrow = f, ncol = 3, 0), c(f, 3, r), dimnames = NULL)

  # a last column of 0s is added to each object
  for (i in 1:r) {
    R[, , i] <- cbind(X[, , i], 0)
  }

  return(R)

}

adjustmentScale <- function(X, r) {
  for (k in 1:r) {
    s <- medland(X[, , k], 2)
    X[, , k] <- X[, , k]/s
  }
  return(X)
}

angurot <- function(R) {
  # The rotation angle associated to rotation matrix R is computed in the range [-pi,pi]
  t = (sum(diag(R)) - 1)/2

  v <- acos(pmin(pmax(t, -1), 1))  # diagonal sum = trace in Matlab

  if (v < (-pi)) {
    v <- (v + (2 * pi))
  } else if (v > pi) {
    v <- (v - (2 * pi))
  }
  return(v)
}

distAllPairsL2 <- function(X) {
  q <- t(X) %*% X
  n <- ncol(q)
  normx <- matlab::repmat(as.matrix(apply(t(X)^2, 1, sum)), 1, n)
  K <- Re(sqrtm(q * (-2) + normx + t(normx)))
  K <- K - (diag(diag(K)))
  return(K)
}

ejerot <- function(R) {

  toler = 1e-16
  for (i in 1:3) {
    for (j in 1:3) {
      if (abs(R[i]) < toler) {
        R[i] <- 0
      }
    }

  }

  # tolerance
  tol = 1e-11
  dif = abs(R - diag(3))
  count = 1

  i = 1
  while (i <= 3) {
    j = 1
    while (j <= 3) {
      if (dif[[i, j]] < tol) {
        count = count + 1
      }
      j = j + 1
    }
    i = i + 1
  }


  if (count < 9) {

    H <- R - diag(3)
    V <- MASS::Null(t(H) %*% H)
  } else {
    V <- matrix(nrow = 1, ncol = 3, 0)
    V[1, 3] <- 1
  }

  for (i in 1:3) {
    if (abs(V[i]) < tol) {
      V[i] <- 0
    }
  }


  if (V[3] < 0) {
    V <- -V
  }

  v <- t(V)
  return(V)
}

escala <- function(X, Y) {
  n <- nrow(X)
  P <- matrix(nrow = n, ncol = n, 1)

  for (i in 1:n) {
    for (j in i + 1:n) {
      if (j <= n) {
        P[i, j] <- norm(t(as.matrix(Y[i, ])) - t(as.matrix(Y[j, ])), "F")/norm(t(as.matrix(X[i,
                                                                                             ])) - t(as.matrix(X[j, ])), "F")
        P[j, i] <- P[i, j]
      }
    }

  }

  return(P)
}

estimorot3D <- function(X, Y, i, j) {
  tol_norm = 1e-11
  dim <- ncol(X)
  A <- matrix(nrow = dim, ncol = dim, 0)
  B <- matrix(nrow = dim, ncol = dim, 0)
  v <- matrix(nrow = 1, ncol = dim, 0)
  w <- matrix(nrow = 1, ncol = dim, 0)
  a3 <- matrix(nrow = 1, ncol = dim, 0)
  b3 <- matrix(nrow = 1, ncol = dim, 0)

  if (i == j) {
    R = diag(3)
  } else {
    A[1, ] = (X[j, ] - X[i, ])/norm(t(as.matrix(X[j, ])) - t(as.matrix(X[i, ])), "F")  #tomo dos landmarks de l mismo especimen
    B[1, ] = (Y[j, ] - Y[i, ])/norm(t(as.matrix(Y[j, ])) - t(as.matrix(Y[i, ])), "F")

    v <- t(as.matrix(vector.cross(A[1, ], X[j, ])))
    w <- t(as.matrix(vector.cross(B[1, ], Y[j, ])))

    normav <- norm(v, "F")
    normaw <- norm(w, "F")

    if ((normav > tol_norm) & (normaw > tol_norm)) {
      A[2, ] <- v/normav
      B[2, ] <- w/normaw
      norm(as.matrix(A[2, ]), "F")
      a3 <- t(as.matrix(vector.cross(A[1, ], A[2, ])))
      A[3, ] <- a3/norm(a3, "F")
      b3 <- t(as.matrix(vector.cross(B[1, ], B[2, ])))
      B[3, ] <- b3/norm(b3, "F")

      R <- t(A) %*% B

    } else if ((normav < tol_norm) & (normaw < tol_norm)) {
      v <- t(as.matrix(vector.cross(A[1, ], B[1, ])))
      w <- v
      normav <- norm(v, "F")

      if (normav < tol_norm) {
        R <- diag(3)
      } else {
        A[2, ] <- v/normav
        B[2, ] <- A[2, ]

        a3 <- t(as.matrix(vector.cross(A[1, ], A[2, ])))
        A[3, ] <- a3/norm(a3, "F")
        b3 <- t(as.matrix(vector.cross(B[1, ], B[2, ])))
        B[3, ] <- b3/norm(b3, "F")


        R <- t(A) %*% B  #we are aligning planes by aligning their orthogonal vectors!!
      }
    } else {
      R <- diag(3)
    }

  }

  return(R)

}

initialFocus <- function(X, f, r) {
  for (k in 1:r) {
    taux <- matrix(nrow = 1, ncol = 3, 0)  #aux. var. to store the actual center
    # taux <- componentwise median of X([,,k]) could be another choice
    taux <- spatialmed_landmark(X[, , k])
    print(taux)
    X[, , k] <- X[, , k] - (matrix(nrow = f, ncol = 1, 1) %*% taux)  # produces object spatial median to be the origin of the coordinate system
  }

  return(X)
}

matrizrot <- function(a) {

  H <- matrix(nrow = 2, ncol = 2, 0)
  H[1, 1] <- cos(a)
  H[1, 2] <- sin(a)
  H[2, 1] <- -sin(a)
  H[2, 2] <- cos(a)

  return(H)
}



matrizrot3D <- function(e, a) {

  AF <- matrix(nrow = 3, ncol = 3, 0)  # to store eigenvectors (as rows)
  Bloque <- matrix(nrow = 3, ncol = 3, 0)  # to store eigenvalues
  R <- matrix(nrow = 3, ncol = 3, 0)

  Bloque[3, 3] <- 1  # the eigenvalue corresponding to the rotation axis
  Bloque[1:2, 1:2] = matrizrot(a)  # the 2D rotation of angle 'a' in the plane orthogonal to axis 'e'
  V = MASS::Null(t(e))
  AF[3, ] <- e  # the axis is placed on the 3rd row
  AF[1:2, ] <- t(V)  #

  R <- t(AF) %*% Bloque %*% AF

  return(R)
}

medianarep <- function(A) {
  ASIND <- sindiagonal(A)  # aux. matrix
  f <- nrow(A)
  aux <- matrix(nrow = f, ncol = 1, 0)

  for (i in 1:f) {
    aux[i, 1] <- median(ASIND[i, ])  #
  }
  m <- median(aux)
  return(m)
}

medland <- function(A, z) {

  f <- nrow(A) #number of landmarks

  b <- matrix(nrow = (f * (f - 1))/2, ncol = 1, 0)
  aux <- matrix(nrow = f, ncol = 1, 0)
  k <- 1
  for (i in 1:f) {
    aux[i, 1] <- norm(t(as.matrix(A[i, ])), "F")  # Frobenius or L2 vector norm
    # print('comienzo con la fila: ') print(i)
    for (j in i + 1:f) {
      # print('comparo con la fila: ') print(j)
      if (j <= f) {
        # print(j)
        b[k, 1] <- norm(t(as.matrix(A[i, ])) - t(as.matrix(A[j, ])), "F")  # distance between pairs of landmarks
        k <- (k + 1)
      }
    }
  }

  if (z == 1) {
    s <- median(aux)
  }
  if (z == 2) {
    s <- median(b)
  }
  if (z == 3) {
    s <- sum(aux)
  }

  return(s)
}

rotation <- function(X, H, Y) {
  eje <- matrix(nrow = 1, ncol = 3, 0)
  tita <- 0

  aux_list <- rot3D(X, Y)
  eje <- aux_list[[1]]
  tita <- aux_list[[2]]

  H <- matrizrot3D(eje, tita)
  X <- X %*% H
  return(list(X, H))
}

scaleSpecimen <- function(X, Y) {
  E <- escala(X, Y)
  ro <- medianarep(E)  # resistant scale factor for the k-th config
  X <- ro * X
  return(X)
}


sindiagonal <- function(B) {

  f <- nrow(B)
  c <- (ncol(B) - 1)
  A <- matrix(nrow = f, ncol = c, 0)

  for (i in 1:f) {
    for (j in 1:c) {
      if (i <= j) {
        A[i, j] <- B[i, j + 1]
      } else {
        A[i, j] <- B[i, j]
      }
    }
  }
  return(A)
}

spatialmed_config <- function(X) {


  # initialization
  n = nrow(X[, , 1])  #number of landmarks
  p = ncol(X[, , 1])  #dimesions
  r = ncol(X[, 1, ])  #number of objects

  A <- matrix(nrow = r, ncol = p, 0)
  M <- matrix(nrow = n, ncol = p, 0)

  w <- matrix(nrow = 1, ncol = r, 1)
  s <- matrix(nrow = 1, ncol = r, 1)

  aux <- matrix(nrow = 1, ncol = p, 0)




  i <- 1
  for (i in 1:n) {
    for (k in 1:r) {
      A[k, ] = X[i, , k]
    }


    aux <- spatialmed_landmark(A)

    M[i, ] <- aux

    aux <- matrix(nrow = 1, ncol = p, 0)


    s <- matrix(nrow = 1, ncol = r, 0)
    A <- matrix(nrow = r, ncol = p, 0)
  }




  return(M)
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

sqrtm <- function(x) {
  z <- x
  for (i in 1:ncol(x)) {
    for (j in 1:nrow(x)) {
      z[i, j] <- sqrt(as.complex(x[i, j]))
    }
  }
  return(z)
}

translation <- function(X, Y, f) {
  TR <- (Y - X)
  t <- matrix(nrow = 1, ncol = 3, 0)
  t <- apply(TR, 2, median)
  X <- X + (matrix(nrow = f, ncol = 1, 1) %*% t(as.matrix(t)))

  return(X)

}

CrossProduct3D <- function(x, y, i = 1:3) {
  # Converts input vectors to three dimensions, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)

  # Indices should be treated cyclically (i.e., index 4 is 'really' index 1, and so on).
  # Index3D() lets us do that, using R's convention of 1-based (rather than 0-based) arrays.
  Index3D <- function(i) (i - 1)%%3 + 1

  # The i-th component of the cross product is: (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat indices cyclically.
  return(x[Index3D(i + 1)] * y[Index3D(i + 2)] - x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

vector.cross <- function(a, b) {

  if (length(a) != 3 || length(b) != 3) {
    a = c(a, 0)
    b = c(b, 0)
  }
  i1 <- c(2, 3, 1)
  i2 <- c(3, 1, 2)
  return(a[i1] * b[i2] - a[i2] * b[i1])

}

rot3D <- function(X, Y) {
  eje <- matrix(nrow = 1, ncol = 3, 0)
  tita <- 0

  f = nrow(X)  # number of landmarks
  c = ncol(X)  # dimensions

  R <- matrix(nrow = c, ncol = c, 0)

  E <- array(matrix(nrow = f, ncol = c, 0), c(f, c, f), dimnames = NULL)


  ejes <- matrix(nrow = f, ncol = c, 0)
  efinal <- matrix(nrow = 1, ncol = c, 0)
  ejaux <- matrix(nrow = 1, ncol = 3, 0)
  Teta <- matrix(nrow = f, ncol = f, 0)
  for (i in 1:f) {
    for (j in 1:f) {
      R <- estimorot3D(X, Y, i, j)
      valEje <- ejerot(R)
      E[j, , i] <- valEje
      E[i, , j] <- valEje
      valAng <- angurot(R)
      Teta[j, i] <- valAng
      Teta[i, j] <- -(valAng) #optimized
    }
  }

  v <- matrix(nrow = 1, ncol = 3, 0)
  h <- 1
  for (i in 1:f) {
    v <- apply(E[, , i], 2, median)
    while ((norm(t(as.matrix(v)), "F") < 0.001) & (h < f)) {
      if (E[h, 3, i] == 1) {
      } else if (E[h, 1, i] < 0) {
        E[h, , i] <- -E[h, , i]
      }
      h <- (h + 1)
      v <- (apply(E[, , i], 2, median))
    }
    h <- 1
    ejes[i, ] <- v/norm(t(as.matrix(v)), "F")
    v <- matrix(nrow = 1, ncol = 3, 0)
  }
  w <- matrix(nrow = 1, ncol = 3, 0)
  for (k in 1:3) {
    w[1, k] <- median(ejes[, k])
  }
  w <- w/norm(w, "F")
  eje <- w
  tita <- medianarep(Teta)
  return(list(eje, tita))
}

