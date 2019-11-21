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

## ----Rho Correction------------------------------------------------------
fileInput <- system.file("testdata", "chen.csv", package="Irescale")
data <- loadFile(fileInput)
rectifiedI<-rectifyIrho(data,1000)
fn = file.path(tempdir(),"output.csv",fsep = .Platform$file.sep)
saveFile(fn,rectifiedI)
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
vI<-resamplingI(distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
statsVI

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
plotHistogramOverlayNormal(vI,statsVI, main=colnames(input$varOfInterest))

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
corrections<-iCorrection(I,vI)
corrections$newI

## ------------------------------------------------------------------------
fileInput <- system.file("testdata", "chen.csv", package="Irescale")
data <- loadFile(fileInput)
distM<-calculateEuclideanDistance(data$data)
vI<-resamplingI(distM,data$varOfInterest,n = 100000)
rectifiedI<- ItoPearsonCorrelation(vI, length(data))
rectifiedI$newI

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
distM<-calculateEuclideanDistance(input$data)
I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
vI<-resamplingI(distM, input$varOfInterest) # This is the permutation
statsVI<-summaryVector(vI)
corrections<-iCorrection(I,vI)
pvalueIscaled<-calculatePvalue(corrections$scaledData,corrections$newI,corrections$summaryScaledD$mean)
pvalueIscaled

## ------------------------------------------------------------------------
library(Irescale)
fileInput<-system.file("testdata", "chen.csv", package="Irescale")
input<-loadFile(fileInput)
resultsChen<-buildStabilityTable(data=input, times=10, samples=10^2, plots=TRUE)

## ------------------------------------------------------------------------
fileInput <- system.file("testdata", "chen.csv", package="Irescale")
data <- loadFile(fileInput)
resultsChen<-buildStabilityTableForCorrelation(data=data,times=10,samples=10^2,plots=TRUE)

