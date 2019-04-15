## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----whole_analysis------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
data<-loadFile(fileInput)
scaledI<-rescaleI(data,samples=1000, scalingUpTo="MaxMin")
fn = file.path(tempdir(),"output.csv",fsep = .Platform$file.sep)
saveFile(fn,scaledI)
if (file.exists(fn)) 
  #Delete file if it exists
  file.remove(fn)

## ------------------------------------------------------------------------
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
head(read.csv(fileInput))

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
head(input$data)
head(input$varOfInterest)

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-"../inst/testdata/chessboard.csv"
input<-loadChessBoard(fileInput)
head(input$data)
head(input$varOfInterest)

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
distM[1:5,1:5]

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-"../inst/testdata/chessboard.csv"
input<-loadChessBoard(fileInput)
distM<-calculateManhattanDistance(input$data)
distM[1:5,1:5]

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
distW<-calculateWeightedDistMatrix(distM)
distW[1:5,1:5]

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
I

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(1000,distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
statsVI

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(1000,distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
plotHistogramOverlayNormal(vI,statsVI, main=colnames(input$varOfInterest))

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(1000,distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
corrections<-iCorrection(I,vI)
corrections$newI

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(1000,distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
corrections<-iCorrection(I,vI)
pvalueIscaled<-calculatePvalue(corrections$scaledData,corrections$newI,corrections$summaryScaledD$mean)
pvalueIscaled

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
resultsChen<-buildStabilityTable(data=input, times=100, samples=1000, plots=TRUE)

