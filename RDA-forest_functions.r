require(vegan)
require(gradientForest)
require(dplyr)
require(stats)
require(utils)

#' Run Gradient Forest (ordination version)
#'
#' Runs gradient forest analysis on an ordination object generated by \code{vegan::capscale()} or \code{vegan::rda()}
#'
#' @param ordination ordination object
#' @param X matrix of predictor variables
#' @param ntrees number of random forest trees to create for each response variable
#' @param keep principal axes of the ordination to analyze. Default is NULL, which means keeping them all.
#' @param ... additional parameters for the gradientForest() function
makeGF=function(ordination,X,ntrees=500,keep=NULL,...) {
  # keep: which PCs to keep in analysis, for example, c(1:25)
  Y=data.frame(scores(ordination,scaling=1,choices=c(1:length(ordination$CA$eig)))$sites)
  # removing constrained axes
  if(length(grep("CAP",colnames(Y)))>0) {Y=Y[,-grep("CAP",colnames(Y))] }
  if(!(is.null(keep))) { Y=Y[,keep] }
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25,...)
  return(gf)
}

#' Run Gradient Forest (simple)
#'
#' This function is a simple wrapper for \code{gradientForest()} function, uses straight-up response matrix \code{Y}.
#'
#' @param Y matrix of response variables
#' @param X matrix of predictor variables
#' @param ntrees number of random forest trees to create during run
#' @param ... additional parameters for the \code{gradientForest()} function
makeGF.simple=function(Y,X,ntrees=500,...) {
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25,...)
  return(gf)
}

#' Predictor selection
#'
#' Selects important predictors based on the fact that they don't decrease in importance at higher \code{mtry} values, unlike non-influential but correlated predictors. Based on simulations by (Strobl et al 2018).
#'
#' @param Y matrix of response variables, or a distance matrix
#' @param X matrix of predictor variables
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param nreps number of spatial bootstrap replicates
#' @param prop.positive.cutoff proportion of replicates in which the predictor must not decrease in raw importance at higher \code{mtry}.
#' @param top.pcs number of top principal components to use for gradientForest
#' @return A list of three items: \code{goodvars} - predictors selected, \code{delta} - change in raw importance of each predictor in each replicate, \code{prop.positive} - proportion of replicates in which each predictor's importance increased at higher \code{mtry}.
mtrySelection=function(Y,X,covariates=NULL,nreps=11,prop.positive.cutoff=0.5,top.pcs=25) {
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  imps=c();delta=list()
  for(i in 1:nreps){
    message("\n\nreplicate ",i)
    rands=rnorm(nrow(Y))
    if(!is.null(covariates)) {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(as.dist(Y)~rands+Condition(as.matrix(covariates)))
      } else {
        ords=rda(Y~rands+Condition(as.matrix(covariates)))
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(as.dist(Y)~rands)
      } else {
        ords=rda(Y~rands)
      }
    }
    message("mtry ",round(ncol(env)/4),"...")
    gf00=makeGF(ords,X,ntrees=1500,mtry=round(ncol(X)/4),keep=c(1:top.pcs))
    message("\nmtry ",round(0.375*ncol(env)),"...")
    gf1=makeGF(ords,X,ntrees=1500,mtry=round(0.375*ncol(X)),keep=c(1:top.pcs))
    ii1=importance(gf1,type="Raw")
    ii0=importance(gf00,type="Raw")[names(ii1)]
    if (i==1) {
      v.order=names(ii1)
    }
    delta[[i]]=ii1[v.order]-ii0[v.order]
  }
  # differences in raw importances
  delta=data.frame(do.call(cbind,delta))
  sdd=stack(delta)
  sdd$var=v.order
  sdd$var=factor(sdd$var,levels=rev(v.order))

  # proportion of positive importance changes
  dc=data.frame(apply(delta,1,function(x){return(sum(x>=0)/length(x))}))
  dc$var=v.order
  dc$var=factor(dc$var,levels=rev(v.order))
  names(dc)[1]="prop.positive"

  goodvars=row.names(dc)[dc$prop.positive>=prop.positive.cutoff]
  return(list(goodvars=goodvars,delta=sdd,prop.positive=dc))
}

#' Spatial Bootstrap
#'
#' Runs gradient forest on \code{nreps} replicates, each time rotating the cloud of datapoints in a random direction and reforming principal axes to be orthogonal to that direction.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param newX new matrix of predictor variables, for transforming into Y-responses based on the fitted model. If omitted, the original \code{X} matrix will be used.
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param nreps number of spatial bootstrap replicates
#' @param top.pcs number of top principal components to use
#' @param mtry \code{mtry} setting. If left unspecified, defaults to \code{N/3} (where \code{N} is the number of predictors in \code{X}.
#' @return A list of three items: \code{turnovers} - predicted Y-responses for \code{newX}, \code{median.importance} - median R2-weighted importance of each predictor across replicates, \code{all.importances} - R2-weighted importance of each predictor in each replicate.
spatialBootstrap=function(Y,X,newX=NULL, covariates=NULL,nreps=25,top.pcs=25,mtry=NULL) {
  #  Y=IBS;X=env[,goodvars];newX=rasters;nreps=3;covariates=covars;top.pcs=25;mtry=NULL
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  if (is.null(newX)) {
    newX=X
  } else {
    newX=newX[,colnames(newX) %in% colnames(X)]
    # making sure there are no values in newX that are outside the range in X
    ranges=apply(X,2,range)
    for (v in colnames(newX)){
      newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
      newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
    }
  }
  imps=c()
  for(i in 1:nreps){
    message("\nreplicate ",i)
    rands=rnorm(nrow(Y))
    if(!is.null(covariates)) {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(as.dist(Y)~rands+Condition(as.matrix(covariates)))
      } else {
        ords=rda(Y~rands+Condition(as.matrix(covariates)))
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(as.dist(Y)~rands)
      } else {
        ords=rda(Y~rands)
      }
    }
    gf1=makeGF(ords,X,ntrees=1500,mtry=mtry,keep=c(1:top.pcs))
    if (i==1) {
      turnovers=predict(gf1,newX)
      v.order=names(importance(gf1))
    } else {
      turnovers=turnovers+predict(gf1,newX)
    }
    imps=c(imps,importance(gf1)[v.order])
  }
  turnovers=turnovers/nreps
  # computing median importances of variables
  dimps=data.frame(imps)
  dimps$variable=names(imps)
  names(dimps)[1]="importance"
  meds=dimps%>%
    group_by(variable)%>%
    summarise(median(importance))
  meds=as.data.frame(meds)
  importances=meds$`median(importance)`
  names(importances)=meds$variable
  # releveling variables according to median importance
  varorder=meds$variable[order(meds[,2],decreasing=T)]
  turnovers=turnovers[,varorder[varorder %in% colnames(turnovers)]]
  importances=importances[varorder]
  dimps$variable=factor(dimps$variable,levels=rev(varorder))
  return(list(turnovers=turnovers,median.importance=importances,all.importances=dimps))
}

