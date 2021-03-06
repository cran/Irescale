% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rectifiedI.R
\name{calculateDistMatrixFromBoard}
\alias{calculateDistMatrixFromBoard}
\title{Calculates the distance in a chessboard-alike structure.}
\usage{
calculateDistMatrixFromBoard(data)
}
\arguments{
\item{data}{is a 2D data structure.}
}
\value{
distM the distance between each cell.
}
\description{
\code{calculateDistMatrixFromBoard} calculates the distance matrix when the field is divided in a matrix shape (rows and columns). This board could have different number of columns for each row.
For example:
\tabular{cccccc}{
1\tab 1\tab 1\tab 1\tab 1\tab 1 \cr
2\tab 2\tab 2\tab 2\tab 2\tab \cr
3\tab 3\tab 3\tab 3\tab 3\tab \cr
4\tab 4\tab 4\tab 4\tab 4\tab \cr
}
The dimension of obtained squared matrix is given by the square of the maximumn dimension of the original matrix. In the previous example, the matrix will have a size of (36,36).
}
\examples{
fileInput <- system.file("testdata", "chessboard.csv", package="Irescale")
data<-loadChessBoard(fileInput)
distM<-calculateEuclideanDistance(data$data)
}
