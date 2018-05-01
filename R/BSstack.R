#' Bootstrap Stacking model builder.
#'
#' Creates a bootstrapped linear stacked set of Random Forest (RF) models given a set of heterogeneous datasets.
#'
#' Required Packages:
#'   dplyr, randomForest, foreach
#'
#' @param T Number of trees for the individual RF models. (int)
#' @param mtry Number of variables available for splitting at each tree node. If a scalar is given then all models use the given values. If a 1D array is given then each individual model uses the given value. If NA then for each model it will be set to Nfeats/3
#' @param nodesize Minimum size of terminal nodes. If a scalar is given then all models use the given values. If a 1D array is given then each individual model uses the given value. By default all models use 5.
#' @param iter The number of time to bootstrap sample the data. (int)
#' @param CV Cross validation to measure mean-absolute error and correlation coefficient, if NA (default) no CV is performed. Otherwise the value gives the number of folds for CV. If CV<2 then leave-one-out CV is performed. CV is performed utilizing the samples that have full record.
#' @param Xn List containing each dataset to be stacked. If not supplied will be generated from X1, X2, ...
#' @param ECHO Bool, enable to provide output to the user in terms of overlapping samples and runtime for CV.
#' @param Y  Nsample x 1 data table of responses for ALL samples. Must have matching rownames with each individual dataset.
#' @param X1 Data table of first dataset to be stacked. Rownames should be contained within Y.
#' @param X2 Data table of second dataset to be stacked. Rownames should be contained within Y.
#' @param ... Further data tables, X3, X4, ..., Xl.
#'
#' @return If CV != null : A list composed of:
#' [1] List containing [1] individual RF models, [2] Nstack +1 weights and [3] feature names for full record samples. This argument is what is used for BSstack_predict
#' [2] Mean-absolute error calculated using cross validation (scalar).
#' [3] Pearson correlation coefficient between actual and predicted values through cross validation (scalar -1<=r<=1).
#' [4] Individual weights calculate for each fold (CV x Nstack+1 matrix).
#' [5] Out of fold predictions for the overlaping samples.
#' [6] Actual values for the overlaping samples.
#' If CV > 1 : Also
#' [7] The fold assignments for the overlapping samples.
#' If CV = null : Only [1] is returned.
#'
#' @export

BSstack <- function(T=50,mtry=NULL,nodesize=5, iter=25, CV=NA,
                    Xn=NULL, ECHO = TRUE,
                    Y, X1, X2, ...) {

  Xl <- list(...)
  # Save datasets to Xn if not given.
  if(is.null(Xn)){
    Nstack <- 2 + length(Xl) # Possible Number of stacked groups
    Xn <- list()

    Xn[[1]] <- X1
    Xn[[2]] <- X2

    if(Nstack > 3){
      for(l in 3:Nstack){
        Xn[[l]] <- Xl[[l-2]]
      }
    }
  } else{
    Nstack <- length(Xn)
  }


  # All feature names
  fNall <- union(colnames(Xn[[1]]),colnames(Xn[[2]]))

  # Overlapping samples
  sNover <- intersect(rownames(Xn[[1]]),rownames(Xn[[2]]))

    # Assign mtry values for each individual RF
    if(is.null(mtry)){
      mtryS <- rep(1,times = Nstack)
      if(is.null(Xn)){
        mtryS[1] <- max(floor(ncol(Xn[[1]])/3),1)
        mtryS[2] <- max(floor(ncol(Xn[[2]])/3),1)
      }
    } else if (length(mtry) == 1){
      mtryS <- rep(mtry,times= Nstack)
    }
    else if (!is.vector(mtry) | length(mtry) != Nstack){
      stop("mtry must be either a scalar or vector with size = number of databases.")
    }
    else{
      mtryS <- mtry
    }

  # Assign nodesize values for each individual RF
  if(is.null(nodesize)){
    nodesizeS <- rep(5,times = Nstack)
  } else if (length(nodesize) == 1){
    nodesizeS <- rep(nodesize,times=Nstack)
  }
  else if (!is.vector(nodesize) | length(nodesize) != Nstack){
    stop("nodesize must be either a scalar or vector with size = number of databases.")
  }
  else{
    nodesizeS <- nodesize
  }

    # Number of samples and feature names for each dataset
    Nl = matrix(0,Nstack,1)
    fNames = list()

    RFn <- list()

    # Assign values for remaining models.
    for (l in 1:Nstack){
      fNall <- union(fNall, colnames(Xn[[l]]))
      sNover <- intersect(sNover, rownames(Xn[[l]]))
      Nl[l] = nrow(Xn[[l]])
      fNames[[l]] = colnames(Xn[[l]])

      if(is.null(mtry)){mtryS[l] <- max(floor(ncol(Xn[[l]])/3),1)}

      RFn[[l]] <- randomForest(Xn[[l]],Y[rownames(Xn[[l]]),1], xtest=NULL, ytest=NULL,ntree=T, mtry = mtryS[l])
    }

  # Find dataset(s) with full record and save them for BS sampling
  nOverl <- 0
  Xfull <- data.frame()

  for (l in 1:Nstack){
    if(setequal(fNall, fNames[[l]])){
      rNames <- rownames(Xfull)
      Xfull <- bind_rows(Xfull,Xn[[l]])

      # Perserve rownames
      if(nOverl > 0){
        rownames(Xfull)[1:nOverl] <- rNames
      }
      rownames(Xfull)[(nOverl+1):(nOverl+Nl[[l]])] <- rownames(Xn[[l]]) # Preserve colnames
      nOverl <- nOverl + Nl[[l]]
    }
  }

  if(length(rownames(Xfull)) == 0 & length(sNover) == 0){
    stop("No valid samples for stacking found.")
  }


  if(length(rownames(Xfull)) > length(sNover)){
    if (ECHO){ print(paste(nrow(Xfull), " Overlapping samples found.")) }
    if(is.na(CV)) {
      w <- BSVerticalStack(T = T, iter = iter, mtry = mtryS, nodesize=nodesizeS, Y = Y, Xfull = Xfull, Xn = Xn, ECHO = FALSE)
      BSmodel <- list("w" = w, "RFn"= RFn, "fNames" = fNames)

      return(list("BSmodel"=BSmodel,"Overlap"=Xfull))
      # Perform LOO cross validation.
    } else if(CV < 2) {
      yPred <- matrix(0,nrow(Xfull),Nstack+1)
      yAct <- matrix(0,nrow(Xfull),Nstack+1)
      BSl <- list()
      wl <- matrix(0,nrow(Xfull),Nstack+1)
      # Go through each overlapping samples and leave it out for testing.
      for(tsti in 1:nrow(Xfull))
      {
        if (ECHO){ print( paste( "Building model for sample: ",tsti)) }

        Xfulltst <- Xfull[tsti,]
        Xfulltrn <- Xfull[-tsti,]
        tstn <- rownames(Xfulltst)

        Xntrn <- Xn
        for (l in 1:Nstack){
          Xntrn[[l]] <- Xn[[l]][setdiff(rownames(Xn[[l]]),tstn),]
        }
        BSl[[tsti]] <- BSstack(T = T, iter = iter, mtry = mtryS, nodesize=nodesizeS, Xn = Xntrn, Y = Y, ECHO = FALSE)
        wl[tsti,] <- BSl[[tsti]]$BSmodel[[1]]
        yPred[tsti,] <- BSstack_predict(BSl[[tsti]]$BSmodel, Xfulltst)
        yAct[tsti,] <- matrix(Y[tstn,1],1,3)
      }

      mae <- colMeans(abs(yPred-yAct))
      r <- diag(cor(yPred, yAct, method = "pearson"))

      BSmodel <- list("w" = colMeans(wl), "RFn"= RFn, "fNames" = fNames)

      return(list("BSmodel"=BSmodel,"Overlap"=Xfull, "mae"=mae,"corr"=r, "wl"=wl, "Pred"=yPred, "Act"=yAct))
    } else if(CV%%1==0 & CV > 2){
      if (CV > nrow(Xfull)){
        stop("ERROR: CV must be < number of overlapping samples.")
      }

      yPred <- matrix(0,nrow(Xfull),Nstack+1)
      yAct <- matrix(0,nrow(Xfull),Nstack+1)
      BSl <- list()
      wl <- matrix(0,CV,Nstack+1)

      # Assign each sample to a fold
      folds <- cut(sample(1:nrow(Xfull)),breaks=CV, labels=FALSE)

      for(n1 in 1:CV){
        if (ECHO){ print( paste( "Building model for fold: ",n1)) }

        tsti <- which(folds==n1,arr.ind=TRUE)
        Xfulltst <- Xfull[tsti,]
        Xfulltrn <- Xfull[-tsti,]
        tstn <- rownames(Xfull)[tsti]

        Xntrn <- Xn
        for (l in 1:Nstack){
          Xntrn[[l]] <- Xn[[l]][setdiff(rownames(Xn[[l]]),tstn),]
        }
        BSl[[n1]] <- BSstack(T = T, iter = iter, mtry = mtryS, nodesize=nodesizeS, Xn = Xntrn, Y = Y, ECHO = FALSE)
        wl[n1,] <- BSl[[n1]]$BSmodel[[1]]
        yPred[tsti,] <- BSstack_predict(BSl[[n1]]$BSmodel, Xfulltst)
        yAct[tsti,] <- matrix(Y[tstn,1],nrow(Xfulltst),3)

      }

      mae <- colMeans(abs(yPred-yAct))
      r <- diag(cor(yPred, yAct, method = "pearson"))


      BSmodel <- list("w"=colMeans(wl), "RFn"=RFn, "fNames"=fNames)

      return(list("BSmodel"=BSmodel,"Overlap"=Xfull, "mae"=mae,"corr"=r, "wl"=wl, "Pred"=yPred, "Act"=yAct, "folds" = folds))
    } else{
      stop("ERROR: CV must be an integer or NA")
    }

    # Horizontal Stacking
  } else{
    if (ECHO){ print(paste(length(sNover), " Overlapping samples found.")) }
    if(is.na(CV)) {
      w <- BSHorizontalStack(T = T, iter = iter, mtry = mtryS, nodesize=nodesizeS, Y = Y, Cf = sNover, Xn = Xn, ECHO = FALSE)
      BSmodel <- list("w" = w, "RFn"= RFn, "fNames"=fNames)

      return(list("BSmodel"= BSmodel, "Overlap"=sNover))
      # Perform LOO cross validation.
    } else if(CV < 2) {
      yPred <- matrix(0,length(sNover),Nstack+1)
      yAct <- matrix(0,length(sNover),Nstack+1)
      BSl <- list()
      wl <- matrix(0,length(sNover),Nstack+1)
      # Go through each overlapping samples and leave it out for testing.
      for(tsti in 1:length(sNover))
      {
        if (ECHO){ print(paste( "Building model for sample: ",tsti)) }

        tstn <- sNover[tsti]
        Xntrn <- Xn
        Xfulltst <- Xn[[1]][tstn,]
        Xntrn[[1]] <- Xn[[1]][setdiff(rownames(Xn[[1]]),tstn),]
        for(n1 in 2:Nstack){
          Xfulltst <- cbind(Xfulltst, Xn[[n1]][tstn,])
          Xntrn[[n1]] <- Xn[[n1]][setdiff(rownames(Xn[[n1]]),tstn), , drop=FALSE]
        }

        BSl[[tsti]] <- BSstack(T = T, iter = iter, mtry = mtryS, nodesize=nodesizeS, Xn = Xntrn, Y = Y, ECHO = FALSE)
        wl[tsti,] <- BSl[[tsti]]$BSmodel[[1]]
        yPred[tsti,] <- BSstack_predict(BSl[[tsti]]$BSmodel, Xfulltst)
        yAct[tsti,] <- matrix(Y[tstn,1],1,3)
      }

      mae <- colMeans(abs(yPred-yAct))
      r <- diag(cor(yPred, yAct, method = "pearson"))

      BSmodel <- list("w" = colMeans(wl),"RFn"= RFn, "fNames"=fNames)

      return(list("BSmodel"=BSmodel,"Overlap"=sNover, "mae"=mae,"corr"=r, "wl"=wl, "Pred"=yPred, "Act"=yAct))
    } else if(CV%%1==0 & CV > 2){
      if (CV > nrow(Xfull)){
        stop("ERROR: CV must be < number of overlapping samples.")
      }

      yPred <- matrix(0,nrow(Xfull),Nstack+1)
      yAct <- matrix(0,nrow(Xfull),Nstack+1)
      BSl <- list()
      wl <- matrix(0,CV,Nstack+1)

      # Assign each sample to a fold
      folds <- cut(sample(1:nrow(Xfull)),breaks=CV, labels=FALSE)

      for(n1 in 1:CV){
        if (ECHO){ print( paste( "Building model for fold: ",n1)) }

        tsti <- which(folds==n1,arr.ind=TRUE)
        Xfulltst <- Xfull[tsti,]
        Xfulltrn <- Xfull[-tsti,]
        tstn <- rownames(Xfull)[tsti]

        Xntrn <- Xn
        for (l in 1:Nstack){
          Xntrn[[l]] <- Xn[[l]][setdiff(rownames(Xn[[l]]),tstn),]
        }
        BSl[[n1]] <- BSstack(T = T, iter = iter, mtry = mtryS, nodesize=nodesizeS, Xn = Xntrn, Y = Y, ECHO = FALSE)
        wl[n1,] <- BSl[[n1]]$BSmodel[[1]]
        yPred[tsti,] <- BSstack_predict(BSl[[n1]]$BSmodel, Xfulltst)
        yAct[tsti,] <- matrix(Y[tstn,1],nrow(Xfulltst),3)

      }

      mae <- colMeans(abs(yPred-yAct))
      r <- diag(cor(yPred, yAct, method = "pearson"))

      BSmodel <- list("w"=colMeans(wl), "RFn"=RFn, "fNames"=fNames)

      return(list("BSmodel"=BSmodel,"Overlap"=sNover, "mae"=mae,"corr"=r, "wl"=wl, "Pred"=yPred, "Act"=yAct, "folds"=folds))
    } else{
      stop("ERROR: CV must be an integer or NA")
    }
  }
}
