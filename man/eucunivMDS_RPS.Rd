% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eucunivMDS_RPS.R
\name{eucunivMDS_RPS}
\alias{eucunivMDS_RPS}
\title{Given a n x n distance matrix D (not necessarily Euclidean) and a initial
 set X0 of n seeds in k dim (that is, an initial n x k matrix), this function finds
 a set of n points in k dimensions X (a final n x k matrix) through a least-squares
 criterion such that the n x n matrix Dk of euclidean distances among these new points X
 is as close as possible to D.}
\usage{
eucunivMDS_RPS(D, k = 2)
}
\arguments{
\item{D}{distance matrix n x n to be approximated}

\item{k}{dimension of output results}
}
\value{
X  A set of n points in k dimensions
}
\description{
Given a n x n distance matrix D (not necessarily Euclidean) and a initial
 set X0 of n seeds in k dim (that is, an initial n x k matrix), this function finds
 a set of n points in k dimensions X (a final n x k matrix) through a least-squares
 criterion such that the n x n matrix Dk of euclidean distances among these new points X
 is as close as possible to D.
}
\examples{
source = array(matrix(nrow = 8,ncol = 3),c(8,3,3),dimnames = NULL)
source[,,1] <- matrix(c(3,0,0,3,0,1,3,1,1,3,1,0,0,0,0,0,0,1,0,1,1,0,1,0)
                   ,nrow = 8,ncol = 3,byrow = TRUE)
source[,,2] <- matrix(c(3, 0 ,0,3, 0, 0.5,3, 1 ,0.75,3 ,1 ,0,0 ,0 ,0,0, 0 ,1,0, 1, 1,0, 1, 0.25)
                     ,nrow = 8,ncol = 3,byrow = TRUE)
source[,,3] <- matrix(c(5, 2 ,1,3, 0, 1.5,3.4, 1 ,1.75,3 ,1 ,0,0 ,0 ,0,0, 2 ,1,0, 3, 1,0, 1, 0.75)
                     ,nrow = 8,ncol = 3,byrow = TRUE)
result <- RPS::robgit_RPS(source, consenso = FALSE)
distance <- RPS::resdistance_RPS(result)
RPS::eucunivMDS_RPS(distance,2)

}
\author{
Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
}
