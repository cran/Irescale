#' @importFrom Rdpack reprompt
#' @importFrom graphics plot points lines hist curve abline text
#' @importFrom grDevices chull topo.colors
#' @importFrom stats dnorm sd quantile median
#' @importFrom utils read.csv write.csv
#' @importFrom sp Polygon
#' @importFrom e1071 skewness kurtosis
#' @importFrom grDevices colorRampPalette



packages <- c("graphics", "rgeos", "sp", "stats", "e1071")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}


library(graphics)
library(rgeos)
library(sp)
library(stats)
library(e1071)
library(testthat)
# library(ggplot2)
# library(tidyverse)



#######################################################
#################LOADING FILES#########################
#######################################################
#'Loads a chessboard or matrix alike input file.
#'
#'\code{loadChessBoard} is used when the input file has a 2D shape, this is a board shape, and it is only one variable of interest.
#'For example:
#'\tabular{cccccc}{
#'1\tab 1\tab 1\tab 1\tab 1\tab 1 \cr
#'2\tab 2\tab 2\tab 2\tab 2\tab \cr
#'3\tab 3\tab 3\tab 3\tab 3\tab \cr
#'4\tab 4\tab 4\tab 4\tab 4\tab \cr
#'}
#'@param fileName the path and file's name to load.
#'@return data frame with two variables, the first variable is a vector with coordinate x (latitude) and y (longitude), the second variable contains the values of the variable of interest.
#'@export loadChessBoard
#'@examples
#' fileInput <- system.file("testdata", "chessboard.csv", package="Irescale")
#' data<-loadChessBoard(fileInput)
loadChessBoard<-function(fileName){
  original<- read.csv(fileName, comment.char = "#", sep = ",", stringsAsFactors = FALSE, na.strings=c("NA","NaN", " ") )

  if(!is.character(original[,1])){
    original[,1]<-as.character(original[,1])
  }
  newRowNames <- make.unique(original[,1])
  rownames(original)<- newRowNames
  original[,1] <- rep(NULL, nrow(original))

  dimOriginal<-dim(original)
  dataS<-list()
  temp1<-unlist(original)
  r<-length(which(!is.na(temp1)))
  print(r)
  data<-matrix(NA,nrow=r,ncol=2)
  varOfInterest<-vector()
  p=1
  for(i in 1:dimOriginal[1]){ #rows
    for(j in 1:dimOriginal[2]){ #columns
      if(!is.null(original[i,j]) & !is.na(original[i,j])){
        data[p,]<-c(i,j)
        varOfInterest[[p]]<-original[i,j]
        p<-p+1
        #varOfInterest<-cbind(varOfInterest,original[i,j])
      }
    }
  }
  dataS[["data"]]<-data
  dataS[["varOfInterest"]]<-varOfInterest
  return(dataS)
}

#'Loads a file with latitude, longitude and variable of interest
#'
#'\code{loadFile} loads the input file with the following format:
#'\itemize{
#'\item Column 1  represents  the sample Id. It has to be Unique.
#'\item Column 2,3 Lat/Long respectively.
#'\item Column 4 and beyond the variables of interest.
#'}
#'
#'@param fileName the file's name and path.
#@examples data <- loadFile(".data/Kiara.csv")
#'@return it returns a data frame with two variables \eqn{data} and \eqn{varOfInterest}. The variable \eqn{data} is a 2D list with the latitude and longitude respectly, while the variable \eqn{varOfInterest} is a matrix with all the variables to calculate and rescale Moran's I.
#'@export loadFile
#'@examples
#'fileInput <- system.file("testdata", "chessboard.csv", package="Irescale")
#' data<-loadFile(fileInput)

loadFile<-function(fileName){
  data <- read.csv(fileName, comment.char = "#", sep = ",", stringsAsFactors = FALSE, na.strings=c("NA","NaN", " ") )
  if(!is.character(data[,1])){
    data[,1]<-as.character(data[,1])
  }
  newRowNames <- make.unique(data[,1])
  rownames(data)<- newRowNames
  data[,1] <- rep(NULL, nrow(data))
  temp1 <- list()
  temp1$data <- data[,1:2]
  rownames(temp1$data)<-rownames(data)
  temp1$varOfInterest <- cbind(data[,3:length(data)])
  colnames(temp1$varOfInterest)<-colnames(data)[3:ncol(data)]
  rownames(temp1$varOfInterest)<-rownames(data)
  return(temp1)
}

##############################################################################################################
#############################CALCULATE DISTANCE###############################################################
##############################################################################################################
#'Transforms a x,y position in a cartesian plane into a position in a 1D array.
#'
#'\code{coor} Transforms a x,y position in a cartesian plane into a position in a 1D array.
#'@param i, the value of the row.
#'@param j, the value of the column.
#'@param size, the maximum between row and columns of the matrix.
#'@return an integer value that represents the position in the array.
#'@export coor
#'@examples
#'pos<-coor(1,1,10)
coor <- function(i,j,size){
  return(i+(j-1)*size)
}

#'Calculates the distance in a chessboard-alike structure.
#'
#'\code{calculateDistMatrixFromBoard} calculates the distance matrix when the field is divided in a matrix shape (rows and columns). This board could have different number of columns for each row.
#'For example:
#'\tabular{cccccc}{
#'1\tab 1\tab 1\tab 1\tab 1\tab 1 \cr
#'2\tab 2\tab 2\tab 2\tab 2\tab \cr
#'3\tab 3\tab 3\tab 3\tab 3\tab \cr
#'4\tab 4\tab 4\tab 4\tab 4\tab \cr
#'}
#'The dimension of obtained squared matrix is given by the square of the maximumn dimension of the original matrix. In the previous example, the matrix will have a size of (36,36).
#'
#'@param data is a 2D data structure.
#'@return distM the distance between each cell.
#'@export calculateDistMatrixFromBoard
#'@examples
#' fileInput <- system.file("testdata", "chessboard.csv", package="Irescale")
#' data<-loadChessBoard(fileInput)
#' distM<-calculateEuclideanDistance(data$data)
calculateDistMatrixFromBoard <-function(data){

  dimData = dim(data)
  m = max(dimData)
  distM = matrix(nrow=m^2, ncol=m^2)
  for(i in 1:dimData[1]){
    for(j in 1:dimData[2]){
      if(!is.na(data[i,j])){
        r <- coor(i,j,m)
        if((j-1)>0){
          if(!is.na(data[i,j-1])){ #left
            c <- coor(i,j-1,m)
            distM[r,c]=1
          }
          if((i-1)>0 && !is.na(data[i-1,j-1])){ #Top,left
            c <- coor(i-1,j-1,m)
            distM[r,c]=1
          }
          if((i+1<=dimData[1]) && !is.na(data[i+1,j-1])){ #Bottom,Left
            c <- coor(i+1,j-1,m)
            distM[r,c]=1
          }
        }
        if((i-1)>0 && !is.na(data[i-1,j])){ #UP
          c <- coor(i-1,j,m)
          distM[r,c]=1
        }
        if((i+1)<=dimData[1] && !is.na(data[i+1,j])){ #Down
          c <- coor(i+1,j,m)
          distM[r,c]=1
        }
        if((j+1)<=dimData[2]){
          if((i-1)>0 && !is.na(data[i-1,j+1])){ #Right,Top
            c <- coor(i-1,j+1,m)
            distM[r,c]=1
          }
          if(!is.na(data[i,j+1])){ #Right
            c <- coor(i,j+1,m)
            distM[r,c]=1
          }
          if((i+1)<=dimData[1] && !is.na(data[i+1,j+1])){ #Right
            c <- coor(i+1,j+1,m)
            distM[r,c]=1
          }
        }
      }
    }
  }

  distM<-distM[, !apply(is.na(distM), 2, all)]
  distM<-distM[!apply(is.na(distM), 1, all),]
  distM[is.na(distM)]=0
  return(distM)
}

#'Given a 2D data structure, it calculates the euclidean distance among all the points.
#'
#'\code{calculateEuclideanDistance} Computes the euclidean distance betwen all pairs of nodes provided in the input vector.
#'
#'Computes the euclidean distance, \eqn{\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2}}, matrix between each pair of points.
#'
#'@param data 2D data structure for latitute and longitute respectively.
#'@return Matrix, of size \eqn{nrow(data) \times nrow(data)}, with the distance between all the pair of points.
#'@export calculateEuclideanDistance
#'@examples
#' fileInput <- system.file("testdata", "chen.csv", package="Irescale")
#' data<-loadFile(fileInput)
#' distM<-calculateEuclideanDistance(data$data)
calculateEuclideanDistance<-function(data){
  if(dim(data)[2]!=2){
    stop("Data Input does not have the requiered dimensions.")
  }
  distM = matrix(data=0,nrow=nrow(data), ncol = nrow(data))
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      distM[i,j]<-sqrt((data[i,1]-data[j,1])^2 + (data[i,2]-data[j,2])^2)
      distM[j,i]<-distM[i,j]
    }
  }
  return(distM)
}

#' Calculates the manhattan distance.
#'
#' \code{calculateManhattanDistance} Calculates the manhattan distance between each pair of nodes.
#' @param data 2D structure with n rows and 2 colums that represents the coordinate in a plane.
#' @return Matrix, of size \eqn{nrow(data) \times nrow(data)}, with the distance between each the pair of points.
#'@export calculateManhattanDistance
#'@examples
#' fileInput <- system.file("testdata", "chessboard.csv", package="Irescale")
#' data<-loadChessBoard(fileInput)
#' distM<-calculateManhattanDistance(data$data)
calculateManhattanDistance<-function(data){
  if(dim(data)[2]!=2){
    stop("Data Input does not have the requiered dimensions.")
  }
  distM = matrix(data=0,nrow=nrow(data), ncol = nrow(data))
  for(i in 1:(nrow(data)-1)){
    for(j in (i+1):nrow(data)){
      distM[i,j]<-abs(data[i,1]-data[j,1]) +abs(data[i,2]-data[j,2])
      distM[j,i]<-distM[i,j]
    }
  }
  return(distM)
}


#'Calculates a weighted representation of the distance matrix.
#'
#'\code{calculateWeightedDistMatrix} The weighted matrix is used as a standardized version of the distance matrix.
#'
#'Computes the similarity matrix of the distance by taking the reciprocal of the distance \eqn{\frac{1}{d}}. A value of Zero is assigned when this value can not be calculated.
#'The whole reciprocal matrix is scaled by dividing each value by the sum of all the elements of the matrix.
#'@param distM, 2D matrix with the distance among all pair of coordinates.
#'@return weighted distance matrix. The sum of this matrix is 1.
#'@export calculateWeightedDistMatrix
#'@examples
#' fileInput <- system.file("testdata", "chen.csv", package="Irescale")
#' data<-loadFile(fileInput)
#' distM<-calculateEuclideanDistance(data$data)
#' distW<-calculateWeightedDistMatrix(distM)
calculateWeightedDistMatrix<-function(distM){
  temp1 <- 1/distM
  temp1[temp1==Inf]<-0
  cSums = sum(temp1)
  return (temp1/cSums)
}

####################################################################################################
#####################################EXTRA ROUTINES ################################################
####################################################################################################

#'Plots the convexhull polygon from the data (latitude, longitude), and calculates the center of the convexhull and its area.
#'
#'\code{convexHull} Computes the area and centroid of the convex hull from the (latitute, longitude) vector.
#'It provides a plot of how the points are dispersed in the field of interest.
#'
#'Consideration for this function:
#'\itemize{
#'\item It makes usage of chull from rgeos and Polygon from graphics.
#'\item The centroid of the polygon is calculated by averaging the vertices of it.
#'\item The shown plot uses the basic \code{plot} command.
#'
#'}
#'@param X dataframe with two colums, latitute and longitude respectively.
#'@param varOfInterest variable of interest to plot. This variable is needed to color the points on the convexhull.
#'@return A vector with two elements, the first element is the area and the second one is the centroid.
#'The centroid is a list of two elements, latitude and longitude that represents the centroid.
#'To have a visual idea of the returned object, it has the following shape \eqn{[area,[latitude,longitude], plotObject]}.
#'@export convexHull
#'@examples
#' fileInput <- system.file("testdata", "chen.csv", package="Irescale")
#' data<-loadFile(fileInput)
#' area_centroid<-convexHull(data$data,data$varOfInterest)
convexHull <- function(X,varOfInterest){
  #todo change latitute, longitude
  #add color label
  #
  shapeX  = dim(X)
  if(shapeX[2]!=2){
    stop("The input vector should have (nX2) shape")
  }
  title="ConvexHull"
  rbPal <- colorRampPalette(c('red','blue'))
  colores<-topo.colors(10)
  if(!missing(varOfInterest)){
    title=colnames(varOfInterest)
    colores = rbPal(10)[as.numeric(cut(varOfInterest,breaks = 10))]
  }

  hpts <- chull(x = X[,1], y = X[,2])
  hpts <- c(hpts, hpts[1])
  xy.coords <- cbind(x = X[,1], y = X[,2])
  chull.coords <- xy.coords[hpts,]
  chull.poly <- Polygon(chull.coords, hole=F)
  chull.area <- chull.poly@area
  chull.centroid <- colMeans(chull.coords)
  #chull.centroid <- colMeans(X)
  # Define the scatterplot
  plot(X, col =colores,main=title, pch=20,cex = 2)
  points(chull.centroid,pch=13, col = "green", cex=3 )
  lines(X[hpts, ])
  return(list(area=chull.area,centroid=chull.centroid))
}




#'p-value calculation.
#'
#'\code{calculatePvalue} calculates a p-value for the null hypothesis.
#'@param sample the vector that will be used as reference.
#'@param value the value of interest.
#'@param mean the mean of interest.
#'@export calculatePvalue
#'@examples
#' fileInput<-system.file("testdata", "chen.csv", package="Irescale")
#' input<-loadFile(fileInput)
#' distM<-calculateEuclideanDistance(input$data)
#' I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
#' vI<-resamplingI(1000,distM, input$varOfInterest) # This is the permutation
#' statsVI<-summaryVector(vI)
#' corrections<-iCorrection(I,vI)
#' pv<-calculatePvalue(corrections$scaledData,corrections$newI,corrections$summaryScaledD$mean)
calculatePvalue<-function(sample,value,mean){
  temp1<-which(sample<=value)
  p<-length(temp1)/length(sample) #number of elements below I
  if(value<mean){
    return(p)
  }
  return(1-p)
}




#'Saves a report with important statistics to describe the sample.
#'
#'\code{saveFile} Saves a csv report with the following columns: Convex Hull Area, Convex Hull Centroid X, Convex Hull Centroid Y, Sample Size, Ichen, Iscaled, pvalue , Mean, MeanScaled, STD DEV, SDScaled, Q_{1\%}, Q_{1\%Scaled}, $Q_{99\%}, Q_{99\%Scaled}, Max, Max_{Scaled}, Min, Min_{Scaled}, Skew,Skew_Scaled, Kutorsis,Kutorsis_Scaled.
#'@param fileName the name  of the file with the path where the CSV file will be saved.
#'@param results is the vector obtained from running the rescaling process over all the variables of interest.
#'@export saveFile
#'@examples
#' fileInput <- system.file("testdata", "chen.csv", package="Irescale")
#' data <- loadFile(fileInput)
#' scaledI<-rescaleI(data,1000)
#' saveFile(paste(tempdir(),"output.csv"),scaledI)
saveFile<-function(fileName,results){
  rList<-data.frame()
  for(i in 1:length(results)){
    row <-c(results[[i]]$convexH$area,
            results[[i]]$convexH$centroid[1],
            results[[i]]$convexH$centroid[2],
            results[[i]]$sampleSize,
            results[[i]]$I,
            results[[i]]$corrections$newI,
            results[[i]]$pvalueI,
            results[[i]]$statsVI$mean,
            results[[i]]$corrections$summaryScaledData$mean,
            results[[i]]$statsVI$median,
            results[[i]]$corrections$summaryScaledData$median,
            results[[i]]$statsVI$sd,
            results[[i]]$corrections$summaryScaledData$sd,
            results[[i]]$statsVI$Q1,
            results[[i]]$corrections$summaryScaledData$Q1,
            results[[i]]$statsVI$Q99,
            results[[i]]$corrections$summaryScaledData$Q99,
            results[[i]]$statsVI$max,
            results[[i]]$corrections$summaryScaledData$max,
            results[[i]]$statsVI$min,
            results[[i]]$corrections$summaryScaledData$min,
            results[[i]]$statsVI$skew,
            results[[i]]$corrections$summaryScaledData$skew,
            results[[i]]$statsVI$kurt,
            results[[i]]$corrections$summaryScaledData$kurt)
    rList<-rbind(rList, row)
  }
  colnames(rList)<-c("Convex Hull Area",
                     "Convex Hull Centroid X",
                     "Convex Hull Centroid Y",
                     "Sample Size",
                     "Ichen",
                     "Iscaled",
                     "P-value",
                     "Mean",
                     "MeanScaled",
                     "Median",
                     "MedianScaled",
                     "Std Dev",
                     "Std Dev Scaled",
                     "Q0.1%",
                     "Q0.1%Scaled",
                     "Q99.9%",
                     "Q99.9%Scaled",
                     "Max",
                     "MaxScaled",
                     "Min",
                     "MinScaled",
                     "Skewness",
                     "SkewnessScale",
                     "Kutorsis",
                     "KutorsisScale")
  rownames(rList)<-names(results)
  write.csv(rList,file=fileName)
}


#' Calculates statistic for the received vector.
#'
#'\code{summaryVector}. Calculates basic statistic of the received vector, like mean, standard deviation, maximum, minimum, 0.1\% and 99.9\% quantile and median.
#'@param vec the vector to calculate the summary.
#'@return a list with mean, standard deviation, maximum, minimum, 0.1\% and 99.9\% quantile and median of the received vector.
#'@export summaryVector
#'@examples
#' inputFileName<-system.file("testdata", "chen.csv", package="Irescale")
#' input<-loadFile(inputFileName)
#' distM<-calculateEuclideanDistance(input$data)
#' I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
#' vI<-resamplingI(1000,distM, input$varOfInterest)
#' statsVI<-summaryVector(vI)
summaryVector<-function(vec){
  q<-quantile(vec,c(.001,.9999))
  return(list(mean=mean(vec),sd=sd(vec),max=max(vec),min=min(vec),Q1=q[1],Q99=q[2],median=median(vec),skew=skewness(vec), kurt = kurtosis(vec)))
}

#' Creates an overlay of the histogram of the data and the theorical normal distribution.
#'
#'\code{plotHistogramOverlayNormal} Overlays the histogram and the theorical normal distribution.
#'@param vec the vector to plot.
#'@param stats the stats obtained from summaryVector.
#'@param density the number of bins for the histogram, The default value is 30.
#'@param main the title of the histogram, The default value is "Histogram".
#'@export plotHistogramOverlayNormal
#'@examples
#' inputFileName<-system.file("testdata", "chen.csv", package="Irescale")
#' input<-loadFile(inputFileName)
#' distM<-calculateEuclideanDistance(input$data)
#' I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
#' vI<-resamplingI(1000,distM, input$varOfInterest)
#' statsVI<-summaryVector(vI)
#' plotHistogramOverlayNormal(vI,statsVI)
plotHistogramOverlayNormal<-function(vec,stats, density=30, main="Histogram"){
  hist(vec,prob=TRUE, col="green",breaks=50,density=density, main=main,xlab = "Moran's I")
  x<-NULL
  curve(dnorm(x=x, mean=stats$median, sd=sd(vec)),col="darkblue", lwd=2, add=TRUE, yaxt="n", xlab = "Moran's I")
  abline(v=stats$median,col="purple",lwd=2)
}


#'Scales a matrix by column.
#'
#'\code{standardizedByColumn} It considers each column independently to scale them.
#'@param M Matrix to be scaled by column.
#'@return a matrix scaled by column.
#'@export standardizedByColumn

standardizedByColumn <- function(M){
  return(apply(M,2,scale))
}

################################################################################################
#######################MORAN'S I Calculation and Corrections####################################
################################################################################################

#' Calculates the Moran's I using the algorithm proposed by Chen \insertCite{chen2009}{Irescale}.
#'
#'\code{calculateMoranI} Moran's I computing method.
#'\deqn{
#' I = varOfInterest^t \times weightedM \times varOfInterest
#'}
#'
#'@param distM the distance matrix. Altough the equation asks for weighted distant matrix, the paramenter that is required is only the distance matrix because this procedure calculate calculates the weighted distance mantrix by itself.
#'@param varOfInterest the variable of interest to calculate Moran's I.
#'@param scaling if the values are previously scaled, set this parameter to False. The default value is TRUE.
#'@return Moran's I
#'@references{
#'\insertAllCited{}
#'}
#'@export calculateMoranI
#' @examples
#' inputFileName<-system.file("testdata", "chen.csv", package="Irescale")
#' input<-loadFile(inputFileName)
#' distM<-calculateEuclideanDistance(input$data)
#' I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
calculateMoranI<-function(distM, varOfInterest, scaling = TRUE){
  #Validation
  dimDistM = dim(distM)
  disVarOfInterest = dim(varOfInterest)

  if(dimDistM[1]!=dimDistM[2]){
    stop("The distance matrix is not squared")
  }
  if(dimDistM[1]!=disVarOfInterest[1]){
    stop("The size of the distance matriz does not match the size of varOfInterest")
  }
  #/Validation
  if(scaling==TRUE){
    varOfInterest<- scale(varOfInterest)
  }
  weightedM = calculateWeightedDistMatrix(distM)
  temp1<-as.matrix(t(varOfInterest))%*%as.matrix(weightedM)%*%as.matrix(varOfInterest)
  return(temp1[1,1])
}

#' Calculates n permutations of the variable of interest to calculate n different I in order to create the \eqn{Null} distribution.
#'
#'\code{resamplingI} Permute n-1 times the values of the variable of interest to calculate a Null distribution for I. It is done n-1, because one order is the original one, to make sure it is included.
#'@param n number of permutations.
#'@param distM the distance matrix. Although the equation requires a weighted distant matrix, the only parameter that will be needed is the distance matrix. This procedure is able to calculate the weighted distance matrix by itself.
#'@param varOfInterest the variable name  or position of the variable we are interested to calculate the spatial autocorrelation.
#'@param scaling The default value is TRUE. However, if the values are previously scaled, this parameter must be set to FALSE..
#'@return A vector with the n calculated Moran's I.
#'@export resamplingI
#'@examples
#' inputFileName<-system.file("testdata", "chen.csv", package="Irescale")
#' input<-loadFile(inputFileName)
#' distM<-calculateEuclideanDistance(input$data)
#' I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
#' vI<-resamplingI(1000,distM, input$varOfInterest)
resamplingI<-function(n,distM,varOfInterest, scaling=TRUE){
  weightedM = calculateWeightedDistMatrix(distM)
  vectorI = NULL
  vectorI = cbind(vectorI,calculateMoranI(distM,varOfInterest, scaling=TRUE))
  for(i in 1:n-1){
    listData<-sample(as.matrix(varOfInterest),replace = FALSE)
    if(scaling==TRUE){
      listData = scale(listData)
    }
    temp1<-as.matrix(t(listData))%*%as.matrix(weightedM)
    vectorI<-cbind(vectorI,(temp1%*%as.matrix(listData))[1,1])
  }
  return(vectorI)
}


#todo
#First: Increase the number of iterations to 100,0000 and instead of using Q1,Q99 use min, max
#Second: Use Q0.001 and Q99.99
#TLC: Run the experiment 100, with 10000 and get an average. \bar(I) <- get the Q0.01, Q99.99
##Get the list of 100 iterations with information, and last part of the report is the average of all that.

#Everything should be run upto 9999 + 1 (The original Chen from population)

#' Scaling process for Moran's I.
#'
#'\code{iCorrection} . consists in centering the I value (I-median) and scaling by the difference between the median and 1st or 99th quantile. The correction is according to the following equation:
#' \deqn{
  #' I = \left\{
  #'\begin{array}{lr}
  #'\frac{(I-median)}{(median - Q1)}& I < median\\
  #'\frac{(I-median)}{(Q99-median)}& I>median\\
  #'\end{array}\right\}
  #' }
#'@param I Moran's I, It could be computed using calculateMoranI function.
#'@param vI the vector obtained by resamplingI.
#'@param statsVI the statistic vector obtained from summaryVector.
#'@param scalingUpTo the rescaling could be done up to the 0.01\% and 99.9\% quantile or max and min values. The two possible options are: "MaxMin", or "Quantile". The default value for this parameter is Quantile.
#'@return rescaled I
#'@export iCorrection
#' @examples
#' inputFileName<-system.file("testdata", "chen.csv", package="Irescale")
#' input<-loadFile(inputFileName)
#' distM<-calculateEuclideanDistance(input$data)
#' I<-calculateMoranI(distM = distM,varOfInterest = input$varOfInterest)
#' vI<-resamplingI(1000,distM, input$varOfInterest)
#' statsVI<-summaryVector(vI)
#' corrections<-iCorrection(I,vI,scalingUpTo="Quantile")
iCorrection<-function(I,vI,statsVI, scalingUpTo="Quantile"){
  if(missing(statsVI)){
    statsVI <- summaryVector(vI)
  }
  upToMin <- statsVI$Q1
  upToMax <- statsVI$Q99
  if(scalingUpTo == "MaxMin"){
    upToMin <- statsVI$min
    upToMax <- statsVI$max
  }
  corrections<-list()
  corrections$scaledData<- vI-statsVI$median
  greaterMeanIndex = which(vI>statsVI$median)
  corrections$scaledData[greaterMeanIndex]<- corrections$scaledData[greaterMeanIndex]/(upToMax-statsVI$median)
  corrections$scaledData[-greaterMeanIndex]<- corrections$scaledData[-greaterMeanIndex]/(statsVI$median-upToMin)
  corrections$summaryScaledData<-summaryVector(corrections$scaledData)
  if(I >= statsVI$median){
    corrections$newI = (I - statsVI$median)/(upToMax-statsVI$median)
  }else{
    corrections$newI = (I - statsVI$median)/ (statsVI$median-upToMin)
  }
  return(corrections)
}




#' Performs the rescale for all the variables in an input file.
#'
#' \code{rescaleI} It executes the whole analysis for all the measurements in the field.
#' \itemize{
#' \item It plots the histogram with the theorical distribution.
#' \item It plots the convexHull for each variable.
#' \item It calcualtes the area and centroid of the convex hull for each variable.
#' \item It calculates the I and rescale it for every variable.
#' \item It returns an object with the computations.
#'}
#'@param data the data frame obtained from \code{loadFile}.
#'@param samples number of permutations for the resampling method.
#'@param scalingUpTo the rescaling could be done up to the 0.01\% and 99.99\% quantile or max and min values. The two possible options are: "MaxMin", or "Quantile". The default value for this parameter is "Quantile"
#'@return An object with I, rescaleI and statistic summary for the inputs without scaling, the same statistics after scaling them, the p-value and the convexhull information
#'@export rescaleI
#'@examples
#' fileInput <- system.file("testdata", "chen.csv", package="Irescale")
#' data <- loadFile(fileInput)
#' scaledI<-rescaleI(data,1000)
rescaleI<-function(data, samples=10000, scalingUpTo ="Quantile"){
results<-list()
for(i in 1:dim(data$varOfInterest)[2]){
  temp1<-colnames(data$varOfInterest)[i]
  record<-list()
  posNA<-which(is.na(data$varOfInterest[,i]))
  if(length(posNA)>0){
    tempData <- data$data[-posNA,]
    tempVarOfInterest <- cbind(data$varOfInterest[-posNA,i])
  }else{
    tempData <- data$data
    tempVarOfInterest <- cbind(data$varOfInterest[,i])
  }
  colnames(tempVarOfInterest)<-c(colnames(data$varOfInterest)[i])
  record$convexH<-convexHull(data$data,tempVarOfInterest)
  distM<-calculateEuclideanDistance(tempData)
  record$sampleSize<-length(tempVarOfInterest)
  record$I<-calculateMoranI(distM,tempVarOfInterest,scaling=TRUE)
  record$vI<-resamplingI(samples,distM, tempVarOfInterest,scaling=TRUE) # This is the permutation
  record$statsVI <- summaryVector(record$vI)
  record$pvalueI<-calculatePvalue(record$vI,record$I,record$statsVI$mean)
  plotHistogramOverlayNormal(record$vI,record$statsVI, main=temp1)
  record$corrections<-iCorrection(record$I,record$vI, scalingUpTo=scalingUpTo)
  plotHistogramOverlayNormal(record$corrections$scaledData,record$corrections$summaryScaledD, main=paste("Scaled",temp1))
  results[[temp1]]<-record
}
return(results)
}

#'Finds how many iterations are necessary to achieve stability in resampling method.
#'
#'\code{buildStabilityTable}  finds how many iterations are necessary to achieve stability in resampling method, plotting in a log scale.
#'@param data data structure after loading the file using \code{loadFile} function
#'@param times the number of times \code{rescaleI} will be executed. The default value is 100.
#'@param samples size of the resampling method. The default value is 1000
#'@param plots to draw the significance plot
#'@return A vector with the average \eqn{\log(samples)} averages I
#'@export buildStabilityTable
buildStabilityTable<-function(data,times = 100, samples=1000, plots=TRUE){
  results<-list()
  for(i in 1:dim(data$varOfInterest)[2]){
    temp1<-colnames(data$varOfInterest)[i]
    record<-list()
    posNA<-which(is.na(data$varOfInterest[,i]))
    if(length(posNA)>0){
      tempData <- data$data[-posNA,]
      tempVarOfInterest <- cbind(data$varOfInterest[-posNA,i])
    }else{
      tempData <- data$data
      tempVarOfInterest <- cbind(data$varOfInterest[,i])
    }
    colnames(tempVarOfInterest)<-c(colnames(data$varOfInterest)[i])
    distM<-calculateEuclideanDistance(tempData)
    record$I<-calculateMoranI(distM,tempVarOfInterest,scaling=TRUE)
    record$corrections<-NULL
    for(n in 1:times){ #up to 100
      record$vI<-resamplingI(samples,distM, tempVarOfInterest,scaling=TRUE) #Ill calculate upto the maximum
      corrections<-NULL
      for(l in 1:log10(samples)){
        record$statsVI <- summaryVector(record$vI[1:10^l]) #this calculation is for 10,10^2....10^n
        corrections<-cbind(corrections,iCorrection(record$I,record$vI[1:10^l])$newI)
      }
      record$corrections<-rbind(record$corrections, corrections) # I only need NewMoransI not everything
    }
    record$corrections<-colMeans(record$corrections)
    if(plots==TRUE){
      plot(record$corrections, main=temp1, type="b", xlab = "log(sample)", ylab="average I", xlim=c(0,log10(samples)+1))
      text(record$corrections, labels=as.character(record$corrections))
    }
    results[[temp1]]<-record$corrections
  }
  return(results)

}

#data<-loadFile("../inst/testdata/chen.csv")
#buildSignificanceTable(data, samples=100000)


# analysisI<-rescaleI("../../DataSets/Kiara/Kiara.csv")
# saveFile("../../DataSets/Kiara/KiaraICalculations.csv",analysisI)


#do a sensitivity analisys
#means I by the number of iterations.
