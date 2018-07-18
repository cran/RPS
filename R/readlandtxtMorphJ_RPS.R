#' Reads a MorphoJ .txt file and returns it
#' as an array of n x k matrices in s dimensions (s=2 or s=3)
#'
#' @param path Path of file
#' @param dim Dimension of the data (2D or 3D).
#'
#' @return A s-dimensional array of n x k matrices and a list of the corresponding object's names
#'
#' @author Guillermo Pacheco, Viviana Ferraggine, Sebastian Torcida
#' @export
readlandtxtMorphJ_RPS <- function(path, dim) {
    result <- read.table(path)  #get data from the file into an array
    s <- nrow(result)  #number of objects
    p <- ncol(result[1, ]) - 1
    p <- p/dim
    M <- array(matrix(nrow = p, ncol = dim, 0), c(p, dim, s), c(2, 1, 3))  #output structure
    names <- list()
    for (i in 1:s) {
        line <- result[i, ]
        name <- as.character(line[[1]])  #gets the names
        names <- c(names, name)
        line <- line[, 2:length(line)]

        m_specimen <- matrix(nrow = p, ncol = dim, 0)
        elem = 1
        for (k in 1:p) {
            coord = 1
            while (coord <= dim) {
                m_specimen[k, coord] <- as.double(line[[elem]])
                elem <- (elem + 1)
                coord = (coord + 1)
            }
        }
        M[, , i] <- m_specimen
    }
    return(list(M, names))
}
