require(dplyr)
require(vegan)
require(gradientForest)
require(stats)
require(utils)
require(cluster)

#' RDAforest-package
#'
#' Random forest analysis of large response matrices, based on ordination scores.
#'
#' @docType package
#'
#' @author Mikhail Matz \email{matz@utexas.edu}
#'
#' @import dplyr
#' @import vegan
#' @import cluster
#' @import gradientForest
#'
#' @name RDAforest-package
NULL


#' Dummify predictors
#'
#' Turns a dataframe containing numerical and categorical predictors into fully numerical
#'
#' @param X matrix of predictor variables
#' @param lastlevel whether to include all levels of a factor, or omit the last one (for lm designs)
#' @return a dataframe with 0,1 values and column names `factor1_level1`,`factor1_level2` etc (for one original column `factor1` with multiple levels). Originally numerical predictors are left unchanged.
#' @export
dummify=function(X,lastlevel=TRUE){
  Xd=list()
  for (ci in 1:ncol(X)) {
    if(is.factor(X[,ci]) | is.integer(X[,ci]) | is.character(X[,ci])) {
      zzzzz=as.factor(X[,ci])
      dum=data.frame(model.matrix(~0+zzzzz))
      if(lastlevel==FALSE) { dum=dum[,-ncol(dum)] }
      colnames(dum)=paste(colnames(X)[ci],colnames(dum),sep="_")
      colnames(dum)=sub("_zzzzz","_",colnames(dum))
      colnames(dum)=sub("^zzzzz","",colnames(dum))
      colnames(dum)=sub("_model.+","",colnames(dum))
      Xd[[ci]]=dum

    }  else {
      Xd[[ci]]=data.frame(X[,ci])
      names(Xd[[ci]])=colnames(X)[ci]
    }
  }
  Xdd=data.frame(do.call(cbind,Xd))
  return(Xdd)
}

#' Sum up importances
#'
#' Sums up importances of original factors that were dummified using `dummify()`
#'
#' @param gf gradient forest result object
#' @param metadata original (non-dummified) metadata
#' @param importance.cutoff do not sum up importances below that level
#' @return a dataframe of importances of original factors in `metadata`.
#' @export
sum_up_importances=function(gf,metadata,importance.cutoff=0){
  ii=data.frame(extendedForest::importance(gf))
  names(ii)="importance"
  ii$var=row.names(ii)
  is=c();ii0=ii[ii$importance>=importance.cutoff,]
  for (f in colnames(metadata)) {
    i2=ii0[grep(paste("^",f,"_",sep=""),ii0$var),]
    if (nrow(i2)==0) {i2=ii0[grep(paste("^",f,sep=""),ii0$var),] }
    is=append(is,sum(i2$importance))
  }
  ii2=data.frame(cbind(importance=is,var=colnames(metadata)))
  ii2$importance=as.numeric(ii2$importance)
  ii2$var=factor(ii2$var,levels=ii2$var[order(ii2$importance)])
  return(ii2)
}


#' Run Gradient Forest (ordination version)
#'
#' This function runs gradient forest analysis on a vegan ordination object. It returns a gradientForest object.
#'
#' @param ordination ordination object made by \code{vegan::capscale()} or \code{vegan::rda()}
#' @param X matrix of predictor variables
#' @param ntrees number of random forest trees to create during run
#' @param keep principal axes of the ordination to analyze. Default is NULL, which means keeping them all.
#' @param ... additional parameters for the gradientForest() function
#' @export
makeGF=function(ordination,X,ntrees=500,keep=NULL,...) {
  require(vegan)
  require(gradientForest)
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
#' @export
makeGF_simple=function(Y,X,ntrees=1500,...) {
  require(vegan)
  require(gradientForest)
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

#' Variable selection (NOTE: this one is outdated! please use mtrySelJack instead)
#'
#' Performs variable selection based on \code{mtry} criterion: variables that are not important by themselves but are correlated with actually important ones diminish in raw importance at higher \code{mtry} setting (Strobl et al 2018). The function runs spatial bootstrap on an ordination and selects variables that do not show decrease in importance at higher \code{mtry} in \code{prop.positive.cutoff} or more replicates.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param nreps number of spatial bootstrap replicates
#' @param ntrees number of trees to construct in each replicate
#' @param prop.positive.cutoff proportion of replicates that must show non-negative importance change at higher \code{mtry} to keep the variable. Default 0.5.
#' @param top.pcs number of top principal components to consider
#' @return A list of three items: \code{goodvars} - variables selected, \code{median.importance} - median importance of each variable, \code{all.importances} - importances inferred across all variables and replicates.
#' @export
mtrySelection=function(Y,X,covariates=NULL,nreps=15,prop.positive.cutoff=0.5,top.pcs=25,ntrees=1500) {
  require(vegan)
  require(dplyr)
  require(gradientForest)
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  imps=c();delta=list()
  for(i in 1:nreps){
    message("\n\nreplicate ",i)
    rands=rnorm(nrow(Y))
    if(!is.null(covariates)) {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(Y~rands+Condition(as.matrix(covariates)))
      } else {
        ords=rda(Y~rands+Condition(as.matrix(covariates)))
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(Y~rands)
      } else {
        ords=rda(Y~rands)
      }
    }
    message("mtry ",round(ncol(X)/4),"...")
    gf00=makeGF(ords,X,ntrees=ntrees,mtry=round(ncol(X)/4),keep=c(1:top.pcs))
    message("\nmtry ",round(0.375*ncol(X))+1,"...")
    gf1=makeGF(ords,X,ntrees=ntrees,mtry=round(0.375*ncol(X))+1,keep=c(1:top.pcs))
    ii1=randomForest::importance(gf1,type=1)
    ii0=randomForest::importance(gf00,type=1)[names(ii1)]
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


#' Variable selection (jackknife version)
#'
#' Performs variable selection based on \code{mtry} criterion: variables that are not important by themselves but are correlated with actually important ones diminish in raw importance at higher \code{mtry} setting (Strobl et al 2018). The function runs ordination jackknife and selects variables that do not show decrease in importance at higher \code{mtry} in \code{prop.positive.cutoff} or more replicates.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param mintry lower \code{mtry} value to use; if omitted, will use 0.25*N of all predictors or 3, whatever is greater.
#' @param maxtry higher \code{mtry} value to use; if omitted, will use 0.375*N+1 of all predictors or 5, whatever is greater.
#' @param oob out-of-bag fraction: the fraction of datapoints to withhold from ordination construction in each replicate
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param nreps number of jackknifing replicates
#' @param ntrees number of trees to construct in each replicate
#' @param prop.positive.cutoff proportion of replicates that must show non-negative importance change at higher \code{mtry} to keep the variable. Default 0.5.
#' @param importance.cutoff lowest R-square importance (at higher mtry) to accept predictor; default 0.01
#' @param top.pcs number of top principal components to consider
#' @return A list of three items: \code{goodvars} - variables selected, \code{median.importance} - median importance of each variable, \code{all.importances} - importances inferred across all variables and replicates.
#' @export
mtrySelJack=function(Y,X,mintry=NULL,maxtry=NULL,nreps=25,oob=0.2,prop.positive.cutoff=0.5,importance.cutoff=0.01,top.pcs=15,covariates=NULL,ntrees=1500) {
  require(vegan)
  require(dplyr)
  require(gradientForest)
  #  Y=sc.fit;X=env.fit;covariates=NULL;nreps=2;oob=0.1;prop.positive.cutoff=0.5;top.pcs=3;ntrees=1500
    if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  imps=c();replicate=c();delta=list()
  for(i in 1:nreps){
    message("\nreplicate ",i)
    ib.keep=sample(1:nrow(Y),round(nrow(Y)*(1-oob)))
    if(!is.null(covariates)) {
      covars.ib=covariates[ib.keep,]
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~Condition(as.matrix(covars.ib)))
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~Condition(as.matrix(covars.ib)))
#        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
        ords=predict(ords0,Y,type='wa',scaling="sites")
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~1)
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
#        ords=predict(ords0,Y,type='wa',scaling="sites")
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~1)
#        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
        ords=predict(ords0,Y,type='wa',scaling="sites")
      }
    }
    if(is.null(mintry)) { mintry=max(3,round(ncol(X)/5)) }
    if(is.null(maxtry)) { maxtry=max(6,round(2*ncol(X)/5)+1) }
    # mintry=min(3,1+round(ncol(X)/4))
    # maxtry=min(6,2+round(0.375*ncol(X)))
    message("mtry ",mintry,"...")
    gf00=makeGF_simple(ords[,c(1:top.pcs)],X,ntrees=ntrees,mtry=mintry)
    message("\nmtry ",maxtry,"...")
    gf1=makeGF_simple(ords[,c(1:top.pcs)],X,ntrees=ntrees,mtry=maxtry)
    ii1=extendedForest::importance(gf1,type="Weighted")
    ii0=extendedForest::importance(gf00,type="Weighted")[names(ii1)]
    if (i==1) {
      v.order=names(ii1)
    }
    delta[[i]]=ii1[v.order]-ii0[v.order]
    imps=c(imps,ii1[v.order])
  }
  # differences in R2-scaled importances
  delta=data.frame(do.call(cbind,delta))
  sdd=stack(delta)
  sdd$var=v.order
  sdd$var=factor(sdd$var,levels=rev(v.order))
  dimps=data.frame(imps)
  dimps$variable=names(imps)
  dimps$variable=factor(dimps$variable,levels=rev(v.order))
  meds=dimps%>%
    group_by(variable)%>%
    summarise(median(imps))
  meds=as.data.frame(meds)
  importances=meds$`median(imps)`
  names(importances)=meds$variable
  # proportion of positive importance changes
  dc=data.frame(apply(delta,1,function(x){return(sum(x>=0)/length(x))}))
  dc$var=v.order
  dc$var=factor(dc$var,levels=rev(v.order))
  names(dc)[1]="prop.positive"
  goodimps=names(importances)[importances>importance.cutoff]
  gooddelta=row.names(dc)[dc$prop.positive>=prop.positive.cutoff]
  goodvars=intersect(goodimps,gooddelta)
  return(list(goodvars=goodvars,delta=sdd,importances=importances,prop.positive=dc))
}

#' Predict multi-column matrix by random forest
#'
#' builds \code{extendedForest::randomForest} model for each column in \code{Y} based on \code{X}
#' then predicts \code{Y} column values based on \code{newX} values of predictors
#'
#' @param Y data frame of columns to predict
#' @param X matrix of predictor variables to build random forest models
#' @param newX new matrix of predictor variable where \code{Y} columns must be predicted.
#' @param extra if newX variable is outside the range of its X values by not more than this fraction of the variable's span, it is set to max or min of X. Values further away are set to NA and these rows are deleted.
#' @param ... other parameters for \code{extendedForest::randomForest}
#' @return list of two items: preds - dataframe of predicted \code{Y} values; goodrows - vector of passing (in-range) rows (T/F, subsetter for plotting)
#' @export
predict_rf=function(Y, X, newX, extra=0, ...) {
  require(extendedForest)
  #  Y=pcs;X=env[,evars];newX=rasters.inrange;extra=0
  # making sure there are no values in newX that are outside the range in X
  newX=newX[,colnames(newX) %in% colnames(X)]
  ranges=apply(X,2,range)
  spans=ranges[2,]-ranges[1,]
  ranges.s=ranges
  ranges.s[2,]=ranges[2,]+spans*extra
  ranges.s[1,]=ranges[1,]-spans*extra
  for (v in colnames(newX)){
    newX[,v][newX[,v]<ranges.s[1,v]]=NA
    newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
    newX[,v][newX[,v]>ranges.s[2,v]]=NA
    newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
  }
  goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
  newX=newX[goodrows,]
  Yp=list()
  for(i in 1:ncol(Y)){
    message("predicting mds ",i)
    pc=Y[,i]
    rf=extendedForest::randomForest(pc~.,data=cbind(pc,X),...)
    Yp[[i]]=predict(rf,newX)
  }
  return(list(preds=data.frame(do.call(cbind,Yp)),goodrows=goodrows))
}

#' Predict multi-column matrix by gradient forest (make turnover curves for clustering)
#'
#' builds \code{extendedForest::randomForest} model for each column in \code{Y} based on \code{X}
#' then predicts \code{Y} column values based on \code{newX} values of predictors
#'
#' @param Y data frame of columns to predict
#' @param X matrix of predictor variables to build random forest models
#' @param newX new matrix of predictor variable where \code{Y} columns must be predicted.
#' @param extra if newX variable is outside the range of its X values by not more than this fraction of the variable's span, it is set to max or min of X. Values further away are set to NA and these rows are deleted.
#' @param ... other parameters for \code{extendedForest::randomForest}
#' @return list of two items: preds - dataframe of predicted \code{Y} values; goodrows - vector of passing (in-range) rows (T/F, subsetter for plotting)
#' @export
predict_gf=function(Y,X,newX,extra=0,...){
  require(gradientForest)
#  Y=pcs;X=env[,evars];newX=rasters
  # making sure there are no values in newX that are outside the range in X
  newX=newX[,colnames(newX) %in% colnames(X)]
  ranges=apply(X,2,range)
  spans=ranges[2,]-ranges[1,]
  ranges.s=ranges
  ranges.s[2,]=ranges[2,]+spans*extra
  ranges.s[1,]=ranges[1,]-spans*extra
  for (v in colnames(newX)){
    newX[,v][newX[,v]<ranges.s[1,v]]=NA
    newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
    newX[,v][newX[,v]>ranges.s[2,v]]=NA
    newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
  }
  goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
  newX=newX[goodrows,]
  gf=makeGF_simple(Y,X)
  pp=predict(gf,newX)
  return(list(preds=pp,goodrows=goodrows))
}


#' Ordination Jackknife
#'
#' Runs gradient forest on \code{nreps} replicates, each time re-building the ordination based on a fraction of the data, projecting the rest of datapoints.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param oob out-of-bag fraction: the fraction of datapoints to withhold from ordination construction in each replicate
#' @param newX new matrix of predictor variable where associated differences in Y must be predicted (averaged across replicates). If omitted, predictions will be made for the original X.
#' @param predict.mode "direct" will use randomForest to predict \code{Y} matrix values. Use this mode to plot geographical maps of local adaptation. "turnover" will predict turnover curves for each variable in \code{X} (i.e., where \code{Y} changes along the \code{X}-variable range). Do NOT use this mode to plot adaptation maps.
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param extra if newX variable is outside the range of its X values by not more than this fraction of the variable's span, it is set to max or min of X. Values further away are set to NA and these rows are deleted.
#' @param nreps number of spatial bootstrap replicates
#' @param top.pcs number of top principal components to consider
#' @param ntrees number of trees to construct in each replicate
#' @param mtry \code{mtry} setting. If left unspecified, defaults to \code{N/3} (where \code{N} is the number of predictors in X)
#' @return A list of three items: \code{predictions} - predicted Y values or turnover curves for \code{newX}, \code{median.importance} - median importance of each variable, \code{all.importances} - all observed importances across all replicates. Table with three columns: \code{importance}, \code{variable}, \code{rep}.
#' @export
ordinationJackknife=function(Y,X,oob=0.2,newX=NULL,predict.mode="direct",covariates=NULL,nreps=25,top.pcs=15,mtry=NULL,ntrees=1500,extra=0) {
#     ntrees=1500;Y=IBS;oob=0.3;X=env[,mm$goodvars];newX=rasters;nreps=3;covariates=covars;top.pcs=5;mtry=NULL
  require(vegan)
  require(dplyr)
  require(gradientForest)
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  if (!is.null(newX)) {
    newX=newX[,colnames(newX) %in% colnames(X)]
    ranges=apply(X,2,range)
    spans=ranges[2,]-ranges[1,]
    ranges.s=ranges
    ranges.s[2,]=ranges[2,]+spans*extra
    ranges.s[1,]=ranges[1,]-spans*extra
    for (v in colnames(newX)){
      newX[,v][newX[,v]<ranges.s[1,v]]=NA
      newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
      newX[,v][newX[,v]>ranges.s[2,v]]=NA
      newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
    }
    goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
    newX=newX[goodrows,]
  } else {
    newX=X
    goodrows=c(1:nrow(X))
    }
  imps=c();replicate=c()

  if(!is.null(covariates)) {
    ords.full=capscale(Y~Condition(as.matrix(covariates)))
  } else {
    ords.full=capscale(Y~1)
  }
  sc.full=scores(ords.full,scaling=1,display="sites",choices=c(1:length(ords.full$CA$eig)))
  for(i in 1:nreps){
    message("\nreplicate ",i)
    ib.keep=sample(1:nrow(Y),round(nrow(Y)*(1-oob)))
    if(!is.null(covariates)) {
      covars.ib=covariates[ib.keep,]
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~Condition(as.matrix(covars.ib)))
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~Condition(as.matrix(covars.ib)))
        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~1)
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~1)
        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
      }
    }
    orr=data.frame(procrustes(sc.full,ords,scale = FALSE)$Yrot)
    gf1=makeGF_simple(orr[,1:top.pcs],X,ntrees=ntrees,mtry=mtry)
    if (i==1) {
          if (predict.mode=="direct") {
            preds=predict_rf(ords[,1:top.pcs],X,newX)$preds
          }
          if (predict.mode=="turnover") {
            preds=predict(gf1,newX)
          }
      v.order=names(extendedForest::importance(gf1,type="Weighted"))
    } else {
        if (predict.mode=="direct") {
          preds=preds+predict_rf(ords[,1:top.pcs],X,newX)$preds
        }
        if (predict.mode=="turnover") {
          preds=preds+predict(gf1,newX)
        }
    }
    imps=c(imps,extendedForest::importance(gf1,type="Weighted")[v.order])
    replicate=c(replicate,rep(i,length(extendedForest::importance(gf1))))
  }
  preds=preds/nreps
  # computing median importances of variables
  dimps=data.frame(imps)
  dimps$replicate=replicate
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
#  turnovers=turnovers[,varorder[varorder %in% colnames(turnovers)]]
  importances=importances[varorder]
  dimps$variable=factor(dimps$variable,levels=rev(varorder))
  return(list(predictions=preds,median.importance=importances,all.importances=dimps,goodrows=goodrows))
}

#' Plot turnover curves returned by \code{ordinationJackknife}
#'
#' @param oj object returned by \code{ordinationJackknife}
#' @param X matrix of predictor variables used to form predictions in \code{oj}
#' @param variable name of the variable to plot
#' @export
plot_turnover=function(oj,X,variable) {
  tt=data.frame(cbind(X,oj$predictions))
  ttt=tt[order(tt[,variable]),]
  plot(ttt[,paste(variable,".1",sep="")]~ttt[,paste(variable)],type="step",xlab="X value",ylab="cumulative Rsquared",main=var2plot,mgp=c(2.3,1,0))
}

#' Plot adaptive neighborhoods
#'
#' Plots a map of differential local adaptation, colored by the first two or three PCs. Can cluster points with similar predicted adaptation, using either the predictions themselves or an additional dataset that is particularly good for clustering (such as predictions from \code{gradientForest}). In the latter case, clusters that contain similar predicted adaptation can be merged.
#' Plots three plots:
#' - hierarchical clustering tree for adaptation clusters with red line for the \code{cluster.merge.level} setting (clusters below this line will be merged);
#' - scatterplot of first two adaptation PCs with \code{envfit} arrows based on the supplied \code{envs} table of environmental variables;
#' - physical map of adaptive neighborhoods, based on supplied spatial coordinates (\code{XY}).
#'
#' @param Y matrix of response variables or a distance matrix
#' @param envs environmental variables to show as \code{envfit} arrows
#' @param xy spatial coordinates of samples in \code{Y}
#' @param nclust number of color clusters to make
#' @param cluster.guide dataset to derive clusters from. Designed to accomodate turnover curves made by \code{predict_gf}
#' @param cluster.merge.level merge clusters that are more similar than this fraction of max cluster distance. Examine the hierarchical clustering tree plotted by this function to set this threshold.
#' @param cluster.labels whether to add cluster numbers to map
#' @param pcs2show how many PCs to derive color scheme from (either 2 or 3). 2 (default) is best for colors to reflect cluster similarity in PCA space.
#' @param coltxt color of arrows and text labels
#' @param flip1 flip flag for PC1 (1 or -1)
#' @param flip2 flip flag for PC2 (1 or -1)
#' @param flip3 flip flag for PC3 (1 or -1)
#' @param scal scaling of envfit arrows: the higher, the shorter the arrows
#' @param jitscale controls distance from arrow point to text label
#' @param rangeExp controls overall X,Y range of the PCA plot (increase to make PCA itself more compact)
#' @param cex text and symbol size control
#' @param ... other \code{plot} options for the final plot (the map)
#' @return (invisible) list of six items: \code{xy} - spatial coordinates; \code{clusters} - clustering indices; \code{PCs} - PC1 and PC2 coordinates; \code{vars} - predictor variable names; \code{arrows} - coordinates of (scaled) arrows on PCA plot, \code{colors} - colors for points.
#' @export
plot_adaptation=function(Y,envs,xy,cluster.guide=NULL,nclust=0,cluster.merge.level=0.333,flip1=1,flip2=1,flip3=1,scal=15,jitscale=0.03,coltxt="black",rangeExp=2,cex=1,cluster.labels=TRUE,pcs2show=2,...){
  require(vegan)
#  require(clv)
#  require(mclust)
#  require(GDAtools)
  require(cluster)
#  Y=rfs;envs=env[rfs0$goodrows,mm$goodvars];xy=XY[rfs0$goodrows,];cluster.guide=pp;flip1=1;flip2=1;flip3=1;cluster.merge.level=0.5;nclust=10;flip2=-1;flip1=1;scal=8;cex=1;jitscale=0.05;rangeExp=1.5;coltxt="black";cluster.labels=TRUE;pcs2show=2
  if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]){ ordination=capscale(Y~1)} else { ordination=rda(Y~1)}
  nc=min(c(length(ordination$CA$eig),3))
  flips=diag(nc)
  diag(flips)=c(flip1,flip2,flip3)[1:nc]
  rd.sc=data.frame(scores(ordination,display="sites",scaling=1,choices=c(1:nc)))
  envft=envfit(ordination,envs,choices=c(1:nc))
  env.arrows=data.frame(scores(envft,display="vectors",scaling=1,choices=c(1:nc)))
  bests=row.names(env.arrows)
  rd.sc=as.matrix(rd.sc) %*% flips
  env.arrows=as.matrix(env.arrows) %*% flips
  pc1 = rd.sc[, 1]
  pc2 = rd.sc[, 2]
  if(pcs2show==2 | nc==2) {
    b <- pc1-pc2
    g <- -pc1
    r <- pc2
  } else {
    pc3=rd.sc[, 3]
    r <- pc1+pc2
    g <- -pc3
    b <- pc3+pc2-pc1
    }
  r <- (r - min(r))/(max(r) - min(r)) * 255
  g <- (g - min(g))/(max(g) - min(g)) * 255
  b <- (b - min(b))/(max(b) - min(b)) * 255
  clusters=NULL
  if(nclust==0) {
    colors=rgb(r,g,b,max=255)
    meds=NULL
    clusters=NULL
  } else {
    if (nclust>0) {
      if(is.null(cluster.guide)) {
        PCs=data.frame(scores(ordination,display="sites",scaling=1,choices=c(1:length(ordination$CA$eig))))
      } else {
        ordi2=rda(cluster.guide~1)
        PCs=data.frame(scores(ordi2,display="sites",scaling=1,choices=c(1:length(ordi2$CA$eig))))
      }
      # message("computing optimal number of clusters...")
      # set.seed(1234)
      # mcc=Mclust(PCs,G=2:nclust)
      # set.seed(NULL)
      # message(mcc$G)
      clPCs = clara(PCs, nclust, sampsize = 500)
      clusters=clPCs$clustering
#      clusters=as.integer(mcc$classification)
#      meds=GDAtools::medoids(dist(PCs),clusters)
      meds=clPCs$i.med
      cdd=cls.scatt.data(Y, clusters)$intercls.centroid
      colnames(cdd)=row.names(cdd)=names(meds)=unique(clusters)
      hc=hclust(as.dist(cdd),method="ave")
      mhi=hc$height[which(hc$height==max(hc$height))]*cluster.merge.level
      if(!is.null(cluster.guide) & cluster.merge.level>0) {
        cl.indices.merged=cutree(hc,h=mhi)
        names(cl.indices.merged)=unique(clusters)
        clusters2=clusters
        clist=sort(unique(clusters))
        names(cl.indices.merged)=clist
        newnames=clist
         seen=c();i=1
        for(i in 1:length(clist)){
           cl=clist[i]
            clusters2[clusters==cl]=cl.indices.merged[cl]
             if(cl.indices.merged[cl] %in% seen) {
               newnames[i]=names(seen)[seen %in% cl.indices.merged[cl]]
             } else {
                 seen=c(seen,cl.indices.merged[cl])
               }
        }
#        names(meds)=cl.indices.merged
        names(meds)=newnames
        clusters=clusters2
        colnames(cdd)=row.names(cdd)=newnames
      }
      hc=hclust(as.dist(cdd),method="ave")
      plot(hc)
      abline(h=mhi,col="red")
    }
    medcolR = c()
    medcolG = c()
    medcolB = c()
    # if(color.by.medoids==TRUE){
    #     medcolR=r[meds]
    #     medcolG=g[meds]
    #     medcolB=b[meds]
    # } else {
        clis=clusters
        for(i in 1:length(unique(clusters))){
          ci = unique(clusters)[i]
          medcolR=c(medcolR,mean(r[clusters==ci]))
          medcolG=c(medcolG,mean(g[clusters==ci]))
          medcolB=c(medcolB,mean(b[clusters==ci]))
          clis[clusters==ci]=i
 #       }
    }
    medcolR <- (medcolR - min(medcolR))/(max(medcolR) - min(medcolR)) * 255
    medcolG <- (medcolG - min(medcolG))/(max(medcolG) - min(medcolG)) * 255
    medcolB <- (medcolB - min(medcolB))/(max(medcolB) - min(medcolB)) * 255
    colors=rgb(medcolR[clis], medcolG[clis], medcolB[clis],max=255)
  }
  xrng <- range(rd.sc[, 1], env.arrows[, 1]/scal) * rangeExp
  yrng <- range(rd.sc[, 2], env.arrows[, 2]/scal) * rangeExp
  jitscale=jitscale*(xrng[2]-xrng[1])
  plot((rd.sc[, 1:2]), xlim = xrng, ylim = yrng, pch = 18, cex = cex, col = colors, asp = 1,axes=FALSE,xlab="",ylab="")
  arrows(rep(0, length(bests)), rep(0,length(bests)), env.arrows[bests,1]/scal, env.arrows[bests,2]/scal, length = 0.0625,col=coltxt)
  jit = rnorm(length(bests),jitscale,jitscale/3)
  text(env.arrows[bests,1]/scal + jit * sign(env.arrows[bests,1]), env.arrows[bests,2]/scal + jit * sign(env.arrows[bests,2]),col=coltxt, labels = bests,cex=cex)
  plot(xy, pch=18,cex = cex, asp = 1, col = colors,...)
  if(length(clusters)>1 & cluster.labels==TRUE) {
    text(xy[meds,], labels = names(meds),col=coltxt,cex=cex)
  }
  invisible(list(xy=xy,clusters=clusters,PCs=rd.sc[,1:2],vars=bests,arrows=env.arrows[bests,1]/scal,colors=colors,medoids=meds))
}

#' Spatial Bootstrap (NOTE: this one is outdated! use ordinationJackknife instead)
#'
#' Runs gradient forest on \code{nreps} replicates, each time rotating the cloud of datapoints in a random direction and reforming principal axes to be orthogonal to that direction.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param newX new matrix of predictor variable where associated differences in Y must be predicted. If omitted, predictions will be made on the same data as the model.
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param nreps number of spatial bootstrap replicates
#' @param top.pcs number of top principal components to consider
#' @param ntrees number of trees to construct in each replicate
#' @param extra fraction of span allowance for range expansion
#' @param mtry \code{mtry} setting. If left unspecified, defaults to \code{N/3} (where \code{N} is the number of predictors in X)
#' @return A list of three items: \code{turnovers} - predicted turnover curves for \code{newX}, \code{median.importance} - median importance of each variable, \code{all.importances} - all observed importances across all replicates. Table with three columns: \code{importance}, \code{variable}, \code{rep}.
#' @export
spatialBootstrap=function(Y,X,newX=NULL, covariates=NULL,nreps=25,top.pcs=25,mtry=NULL,ntrees=1500,extra=0) {
  require(vegan)
  require(dplyr)
  require(gradientForest)
  #  Y=IBS;X=env[,goodvars];newX=rasters;nreps=3;cov
  if (is.null(newX)) {
    newX=X
  } else {
    newX=newX[,colnames(newX) %in% colnames(X)]
    # making sure there are no values in newX that are outside the range in X
    ranges=apply(X,2,range)
    spans=ranges[2,]-ranges[1,]
    ranges.s=ranges
    ranges.s[2,]=ranges[2,]+spans*extra
    ranges.s[1,]=ranges[1,]-spans*extra
    for (v in colnames(newX)){
      newX[,v][newX[,v]<ranges.s[1,v]]=NA
      newX[,v][newX[,v]>ranges.s[2,v]]=NA
      newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
      newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
    }
    goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
    newX=newX[goodrows,]
  }
  imps=c();replicate=c()
  for(i in 1:nreps){
    message("\nreplicate ",i)
    rands=rnorm(nrow(Y))
    if(!is.null(covariates)) {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(Y~rands+Condition(as.matrix(covariates)))
      } else {
        ords=rda(Y~rands+Condition(as.matrix(covariates)))
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2]) {
        ords=capscale(Y~rands)
      } else {
        ords=rda(Y~rands)
      }
    }
    gf1=makeGF(ords,X,ntrees=ntrees,mtry=mtry,keep=c(1:top.pcs))
    if (i==1) {
      turnovers=predict(gf1,newX)
      v.order=names(extendedForest::importance(gf1))
    } else {
      turnovers=turnovers+predict(gf1,newX)
    }
    imps=c(imps,extendedForest::importance(gf1)[v.order])
    replicate=c(replicate,rep(i,length(extendedForest::importance(gf1))))
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

