% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rectifiedI.R
\name{transformImageToList}
\alias{transformImageToList}
\title{Transforms the image in the object need it to run the analysis.}
\usage{
transformImageToList(im)
}
\arguments{
\item{im}{cimg object.}
}
\description{
\code{transformImageToList} transforms the image into a list with two variables, data and varOfInterest, which are the identificators needed to execute the rectification.
}
\examples{
fileInput <- system.file("testdata", "imageGray.png", package="Irescale")
img<-loadSatelliteImage(fileInput)
data<-transformImageToList(img)
}
