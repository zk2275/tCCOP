## =====================================================================
## tCOPP_functions.R  -- shared function library for the tCOPP simulations
## ---------------------------------------------------------------------
## Source this ONCE per session, then run any simulation file's
## simulation/plotting part:   source("tCOPP_functions.R")
##
## Contains the 35 functions that were byte-identical across all five
## simulation .Rmd files, PLUS the two MuSCI_ConformalCI* functions
## (SuperLearner version, used by 4 of 5 files). The B==1 path has been
## fixed (see "B==1 FIX" markers) so B genuinely controls the number of
## splits for ANY method:  B>1 -> majority-vote aggregation across splits,
## B==1 -> a single split.
##
## NOT included (deliberately, because they legitimately differ per file):
##   - genY_singlestage_ex : the data-generating process is figure-specific
##                           and must stay in each simulation file.
## NOTE on MuSCI: tCOPP_Fig1 uses an "[Optimized]" glm rewrite that is
## numerically identical to the SuperLearner version here whenever
## m.Library = ps.Library = "SL.glm" (true at every call site). Keep that
## file's own optimized MuSCI if you want its speed; otherwise this shared
## copy is equivalent.
## =====================================================================

library(MASS)
library(numDeriv)
library(rootSolve)
library(osqp)
library(densratio)
library(dplyr)
library(glmnet)
library(quantreg)
library(randomForest)
library(grf)
library(SuperLearner)

MuSCI_ConformalCI_MuS <- function(X,Y,X_ex,Y_ex,T,T_ex,Deta,Xtest,method,B,trainprop,psfun,psparams,outfun,outparams,quantiles,type,side,weighted,
####                  conf.score="ASR",
                  conf.score,
                  ps.Library="SL.glm",
                  m.Library="SL.glm",
                  p1,#0.5
                  seed=10012) {
  YBleft <- NULL
  YBright <- NULL
  Nts <- NULL
  r0<-NULL
  r1<-NULL
  X_copy = data.frame(X)
  Y_copy = Y
  X_ex_copy = X_ex
  Y_ex_copy = Y_ex
  X_copy$Site=1

#random sample from binary varibale(A) B times
  for (b in 1:B) {
  n <- length(Y)
  #split samples for training and calibration
  YY <- Y
  YY[!T] <- NA
  T <- as.numeric(!is.na(YY))
  res <- SampleSplit(YY,trainprop)
  #res$trainid: index for all training
  #res$trainid1: index for all training with T=1
  #!!!need to select T=1 first and then select trainid1
  #res$trainid0: index for all training with T=0

  #estimate the propensity score
  Xtrain <- X_copy[res$trainid, ,drop=FALSE]
  Ttrain <- T[res$trainid]
  PSFit_T <- PropScore(Ttrain,Xtrain[,1:(ncol(Xtrain)-1)],psfun, psparams)

  #estimate weights for testing data
  ps_T_test <- PSFit_T(Xtest)
  Xtest <- as.matrix(Xtest)
  ps_D_test <- exp(cbind(1,Xtest)%*%Deta)/(1+exp(cbind(1,Xtest)%*%Deta))
  ps_test <- ps_T_test*(1-ps_T_test)/(ps_T_test*(1-ps_D_test)+ps_D_test*(1-ps_T_test))


  #estimate the propensity score for binary A
  tp_T <- PSFit_T(X_copy[,1:(ncol(X_copy)-1)])
  tp_D <- exp(as.matrix(cbind(1,X_copy[,1:(ncol(X_copy)-1)]))%*%Deta)/(1+exp(as.matrix(cbind(1,X_copy[,1:(ncol(X_copy)-1)]))%*%Deta))
  tp_A <- tp_D*(1-tp_T)/(tp_D*(1-tp_T)+tp_T*(1-tp_D))


    A <- rbinom(n = n, size = 1, prob=tp_A)
    ##print("Len:")
    ##print(length(T)[1])
    ##print(n)
    ##print(length(A))
    O <- I(T==A)

    Nts[b] <- sum(O)
    #quantile regression fitting after selecting A=T
    X1 <- X_copy[res$trainid, ,drop=FALSE]
    ####X1 <- X1[,1:(ncol(X1)-1)]
    Y1 <- Y_copy[res$trainid]
    #select O=1 from training data
    ind <- which(O[res$trainid] == 1)
    Xtrain <- X1[ind, ,drop=FALSE]
    Ytrain <- Y1[ind]
    ####OutputFit <- QuantileReg(Ytrain,Xtrain,outfun,outparams,quantiles)
    ####Yhat_test <- OutputFit(Xtest)
    ###sample and fit for train
    ####sample for val
    if(method == "Sweight" | method == "BootSweight"){
      #data for validation
      X0 <- X_copy[-res$trainid, ,drop=FALSE]
      Y0 <- Y_copy[-res$trainid]
      ind <- which(O[-res$trainid] == 1)
      Xval <- X0[ind, ,drop=FALSE]
      Yval <- Y0[ind]
      ####ps_T_val <- PSFit_T(as.data.frame(Xval))
      ####Yhat_val <- OutputFit(Xval)

      ####Xval <- as.matrix(Xval)
      ####ps_D_val <- exp(cbind(1,Xval)%*%Deta)/(1+exp(cbind(1,Xval)%*%Deta))
      ####ps_val <- ps_T_val*(1-ps_T_val)/(ps_T_val*(1-ps_D_val)+ps_D_val*(1-ps_T_val))

      ####subsample wttl, wttl is for both train and val
      ####wttl0 = wttl[-res$trainid]
      ####wttl_val = wttl0[ind]
    }
    #else
    s=max(X_ex$Site)
    for (k in 2:max(X_ex$Site)){
      Site_k_idx = which(X_ex$Site==k)
      X_temp=X_ex[Site_k_idx,]
      T_temp=T_ex[Site_k_idx]
      Y_temp=Y_ex[Site_k_idx]
      n <- length(Site_k_idx)
  #split samples for training and calibration
      YY <- Y_temp
      YY[!T_temp] <- NA
      T_temp <- as.numeric(!is.na(YY))
      res <- SampleSplit(YY,trainprop)
  #res$trainid: index for all training
  #res$trainid1: index for all training with T=1
  #!!!need to select T=1 first and then select trainid1
  #res$trainid0: index for all training with T=0

  #estimate the propensity score
      Xtrain_p <- X_temp[res$trainid, ,drop=FALSE]
      Ttrain_p <- T_temp[res$trainid]
      PSFit_T <- PropScore(Ttrain_p,Xtrain_p[,1:(ncol(Xtrain_p)-1)],psfun, psparams)

  #estimate weights for testing data
      ps_T_test <- PSFit_T(Xtest)
      ####Xtest <- as.matrix(Xtest)
      ps_D_test <- exp(as.matrix(cbind(1,Xtest))%*%Deta)/(1+exp(as.matrix(cbind(1,Xtest))%*%Deta))
      ps_test <- ps_T_test*(1-ps_T_test)/(ps_T_test*(1-ps_D_test)+ps_D_test*(1-ps_T_test))


  #estimate the propensity score for binary A
      tp_T <- PSFit_T(X_temp[,1:(ncol(X_temp)-1)])
      tp_D <- exp(as.matrix(cbind(1,X_temp[,1:(ncol(X_temp)-1)]))%*%Deta)/(1+exp(as.matrix(cbind(1,X_temp[,1:(ncol(X_temp)-1)]))%*%Deta))
      tp_A <- tp_D*(1-tp_T)/(tp_D*(1-tp_T)+tp_T*(1-tp_D))


      A <- rbinom(n = n, size = 1, prob=tp_A)
    ##print("Len:")
    ##print(length(T)[1])
    ##print(n)
    ##print(length(A))
      O <- I(T_temp==A)

      Nts[b] <- sum(O)
    #quantile regression fitting after selecting A=T
      X1 <- X_temp[res$trainid, ,drop=FALSE]
      Y1 <- Y_temp[res$trainid]
    #select O=1 from training data
      ind <- which(O[res$trainid] == 1)
      Xtrain <- rbind(Xtrain,X1[ind, ,drop=FALSE])
      Ytrain <- c(Ytrain,Y1[ind])
    ####OutputFit <- QuantileReg(Ytrain,Xtrain,outfun,outparams,quantiles)
    ####Yhat_test <- OutputFit(Xtest)
    ###sample and fit for train
      ####sample for val
    if(method == "Sweight" | method == "BootSweight"){
      #data for validation
      X0 <- X_temp[-res$trainid, ,drop=FALSE]
      Y0 <- Y_temp[-res$trainid]
      ind <- which(O[-res$trainid] == 1)
      Xval <- rbind(Xval,X0[ind, ,drop=FALSE])####20240911
      Yval <- c(Yval,Y0[ind])
      ####ps_T_val <- PSFit_T(as.data.frame(Xval))
      ####Yhat_val <- OutputFit(Xval)

      ####Xval <- as.matrix(Xval)
      ####ps_D_val <- exp(cbind(1,Xval)%*%Deta)/(1+exp(cbind(1,Xval)%*%Deta))
      ####ps_val <- ps_T_val*(1-ps_T_val)/(ps_T_val*(1-ps_D_val)+ps_D_val*(1-ps_T_val))

      ####subsample wttl, wttl is for both train and val
      ####wttl0 = wttl[-res$trainid]
      ####wttl_val = wttl0[ind]
    }
    #else
    }####end for for k in 1:s

  ###############################################################################
  ##################            Data pre-processing           ###################
  ###############################################################################
  ####colnames(data.train)[colnames(data.train)==T.name] <- "Site"
  ####site.names <- unique(data.train$Site) # extract names of sites
  ####s <- length(site.names) # number of sites
  ####s=1
  ####non.tgt.site.names <- site.names[site.names%in%tgt.site.name==F] # names of non-target sites

  ####data.train$Site[data.train$Site==tgt.site.name] <- 999
  ####Xtrain$Site <- 999
  ###for(k in 2:s) {
  ###  data.train$Site[data.train$Site==non.tgt.site.names[k-1]] <- k*999
  ###}
  ####Xtrain$Site <- Xtrain$Site/999
  ####covar.names <- names(Xtrain)[names(Xtrain)%in%c(Y.name, "Site")==F]
  D1 <- D2 <- data.frame()
  set.seed(seed)

  D11=data.frame(Xtrain)
  D21=data.frame(Xval)
  Xtest2=data.frame(Xtest)
  Xtest2$"Site"=1
  colnames(D11)=colnames(D21)=colnames(Xtest2)=c(paste0("cn",1:(ncol(D11)-1)),"Site")####20240911
  covar.names=colnames(D21)[1:(ncol(D21)-1)]
  ###############################################################################
  ##################            Pooled sample method           ##################
  ###############################################################################
  ####D1$R <- as.numeric(is.na(D1[,Y.name])==F) # define the missingness indicator
  ####D2$R <- as.numeric(is.na(D2[,Y.name])==F)
  ####D11 <- D1[D1$R==1,]
  ####D21 <- D2[D2$R==1,]

  ####y = D1$R
  ####x = D1[,covar.names]

  y=c(rep(1,nrow(rbind(D11,D21))),rep(0,nrow(Xtest2)))
  x=rbind(D11,D21,Xtest2)
  x=x[,covar.names]
  ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
  pred = as.numeric(predict(ps.fit, D21[,covar.names])$library.predict)
  prob = pmax(pmin(pred, 0.9), 0.1)
  eta = (1-prob)/prob

  D11all=cbind(D11,Ytrain)
    D21all=cbind(D21,Yval)
    colnames(D11all)[ncol(D11all)]="Yname"
    colnames(D21all)[ncol(D21all)]="Yname"
    Y.name="Yname"
    D11=D11all
    D21=D21all
    res2 <- SampleSplit(Xtest2[,1],trainprop)
    D10=data.frame(cbind(Xtest2[res2$trainid,],rep(NA,nrow(Xtest2[res2$trainid,]))))####20240911
    colnames(D10)=colnames(D11)
    D1=rbind(D11,D10)
    D20=data.frame(cbind(Xtest2[-res2$trainid,],rep(NA,nrow(Xtest2[-res2$trainid,]))))####20240911
    colnames(D20)=colnames(D21)
    D2=rbind(D21,D20)
    colnames(D1)=colnames(D11)
    colnames(D2)=colnames(D21)
    D1$R=D2$R=0
    D1$R[1:nrow(D11)]=1
    D2$R[1:nrow(D21)]=1
    covar.names=colnames(D21)[1:(ncol(D21)-2)]####20240911

  if(conf.score=="ASR") {
    ##print("Y")
    ##print(head(D11[,Y.name]))
    CFS <- ASR(W=data.matrix(D11[,covar.names]),
               Ytemp=D11[,Y.name],
               W.pred=data.matrix(D21[,covar.names]),
               Y.pred=D21[,Y.name])
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  ##print("flag1")
  if(conf.score=="localASR") {
    CFS <- localASR(W=data.matrix(D11[,covar.names]),
                    Ytemp=D11[,Y.name],
                    W.pred=data.matrix(D21[,covar.names]),
                    Y.pred=D21[,Y.name])
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    CFS <- CQR(data=D11, Ytemp=D11[,Y.name],
               data.pred=D21, Y.pred=D21[,Y.name],
               rq.form=rq.form, alpha=alpha)
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }

  # [Optimized] Direct glm + cached m(theta)
  m <- function(theta) {
    y <- as.numeric(CFS1 <= theta)
    x <- D11[, covar.names]
    fit <- glm(y ~ ., data = data.frame(y = y, x), family = binomial())
    p_CFS <- as.numeric(predict(fit, newdata = D2[, covar.names], type = "response"))
    return(p_CFS)
  }
  f <- function(theta){
    m_val <- m(theta)
    sum(m_val[D2$R==0]-(1-alpha))/nrow(D2[D2$R==0,]) +
       sum(eta*(I(CFS2<=theta)-m_val[D2$R==1]))/nrow(D2[D2$R==1,])
  }
  r_0 <- c(CFS1, CFS2)
  rhat_pool = try(uniroot(f,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
  if(class(rhat_pool)!= "try-error"){
    rhat_pool = rhat_pool
  } else if (f(max(r_0)/2)>0) {
    rhat_pool = try(uniroot(f,c(0,max(r_0)/2),tol=1e-5)$root)
  } else {
    rhat_pool = try(uniroot(f,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
  }
  ###############################################################################
  ##################         Use target site data only         ##################
  ###############################################################################
  ####D1_s1 <- D1[D1$Site==1, ]
  ####D2_s1 <- D2[D2$Site==1, ]
  ####D11_s1 <- D1_s1[D1_s1$R==1, ]
  ####D21_s1 <- D2_s1[D2_s1$R==1, ]
  D1_s1 <- D1[D1$Site==1, ]
  D2_s1 <- D2[D2$Site==1, ]
  D11_s1 <- D1_s1[D1_s1$R==1, ]
  D21_s1 <- D2_s1[D2_s1$R==1, ]
  #D1_s1 <- D1
  #D2_s1 <- D2
  #D11_s1 <- D11
  #D21_s1 <- D21

  y = D1_s1$R
  x = D1_s1[,covar.names]
  ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
  pred = as.numeric(predict(ps.fit, D21_s1[,covar.names])$library.predict)
  prob = pmax(pmin(pred, 0.9), 0.1)
  eta_s1 = (1-prob)/prob

  if(conf.score=="ASR") {
    CFS <- ASR(W=data.matrix(D11_s1[,covar.names]),
               Ytemp=D11_s1[,Y.name],
               W.pred=data.matrix(D21_s1[,covar.names]),
               Y.pred=D21_s1[,Y.name])
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  if(conf.score=="localASR") {
    CFS <- localASR(W=data.matrix(D11_s1[,covar.names]),
                    Ytemp=D11_s1[,Y.name],
                    W.pred=data.matrix(D21_s1[,covar.names]),
                    Y.pred=D21_s1[,Y.name])
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    CFS <- CQR(data=D11_s1, Ytemp=D11_s1[,Y.name],
               data.pred=D21_s1, Y.pred=D21_s1[,Y.name],
               rq.form=rq.form, alpha=alpha)
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }

  # [Optimized] Direct glm + cached m_s1
  m_s1 <- function(theta) {
    y <- as.numeric(CFS1_s1 <= theta)
    x <- D11_s1[, covar.names]
    fit <- glm(y ~ ., data = data.frame(y = y, x), family = binomial())
    p_CFS <- as.numeric(predict(fit, newdata = D2_s1[, covar.names], type = "response"))
    return(p_CFS)
  }

  f_s1 <- function(theta){
    m_s1_val <- m_s1(theta)
    sum(m_s1_val[D2_s1$R==0]-(1-alpha))/nrow(D2_s1[D2_s1$R==0,])  +
       sum(eta_s1*(I(CFS2_s1<=theta)-m_s1_val[D2_s1$R==1]))/nrow(D2_s1[D2_s1$R==1,])
  }
  r_0 <- c(CFS1_s1, CFS2_s1)
  rhat_tgSt = try(uniroot(f_s1,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
  if(class(rhat_tgSt)!= "try-error"){
    rhat_tgSt = rhat_tgSt
  } else if (f_s1(max(r_0)/2)>0) {
    rhat_tgSt = try(uniroot(f_s1,c(0,max(r_0)/2),tol=1e-5)$root)
  } else {
    rhat_tgSt = try(uniroot(f_s1,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
  }
  ####
  r0=c(r0,rhat_tgSt)
print("flag")
print("r0")
print(r0)
print("s")
print(s)
###############################################################################
  ##############       Federated weights and equal weights        ###############
  ###############################################################################
  chi <- rhat_ntgt <- c()
  omega <- eta <- phi <- CFS1 <- CFS2 <- list()
  D10_s1 <- D1_s1[D1_s1$R==0, ]
  D20_s1 <- D2_s1[D2_s1$R==0, ]
  eta=list()
  for(k in 2:s) {
    D1_sk <- D1[D1$Site==k, ]
    D2_sk <- D2[D2$Site==k, ]
    D11_sk <- D11[D11$Site==k, ]
    D21_sk <- D21[D21$Site==k, ]

    y = D1_sk$R
    x = D1_sk[,covar.names]
    ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
    pred = as.numeric(predict(ps.fit,  D21_sk[,covar.names])$library.predict)
    prob = pmax(pmin(pred, 0.9), 0.1)
    eta_sk = (1-prob)/prob
    ####print("eta_sk")
    ####print(eta_sk)
    eta[[k-1]] <- eta_sk

    if(conf.score=="ASR") {
      CFS <- ASR(W=data.matrix(D11_sk[,covar.names]),
                 Ytemp=D11_sk[,Y.name],
                 W.pred=data.matrix(D21_sk[,covar.names]),
                 Y.pred=D21_sk[,Y.name])
      CFS1_sk <- CFS[[1]]; CFS2_sk <- CFS[[2]]
    }
    if(conf.score=="localASR") {
      CFS <- localASR(W=data.matrix(D11_sk[,covar.names]),
                      Ytemp=D11_sk[,Y.name],
                      W.pred=data.matrix(D21_sk[,covar.names]),
                      Y.pred=D21_sk[,Y.name])
      CFS1_sk <- CFS[[1]]; CFS2_sk <- CFS[[2]]
    }
    if(conf.score=="CQR") {
      rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
      CFS <- CQR(data=D11_sk, Ytemp=D11_sk[,Y.name],
                 data.pred=D21_sk, Y.pred=D21_sk[,Y.name],
                 rq.form=rq.form, alpha=alpha)
      CFS1_sk <- CFS[[1]]; CFS2_sk <- CFS[[2]]
    }
    CFS1[[k-1]] <- CFS1_sk
    CFS2[[k-1]] <- CFS2_sk

    D10_s1M <- as.matrix(D10_s1)
    D11_skM <- as.matrix(D11_sk)
    D20_s1M <- as.matrix(D20_s1)
    D21_skM <- as.matrix(D21_sk)
    ratio <- Estimate_omega_np(x=D11_skM[,covar.names],
                               x_target=D10_s1M[,covar.names],
                               x.pred=D21_skM[,covar.names],
                               x_target.pred=D20_s1M[,covar.names]) %>%
      as.numeric()
    omega[[k-1]] <- ratio

    # [Optimized] Direct glm + cached m_s1 and m_sk
    m_sk <- function(theta) {
      y <- as.numeric(CFS1_sk <= theta)
      x <- D11_sk[, covar.names]
      fit <- glm(y ~ ., data = data.frame(y = y, x), family = binomial())
      p_CFS_sk <- as.numeric(predict(fit, newdata = D2_sk[, covar.names], type = "response"))
      return(p_CFS_sk)
    }

    f_sk = function(theta){
      m_s1_val <- m_s1(theta)
      m_sk_val <- m_sk(theta)
      sum((m_s1_val[D2_s1$R==0]-(1-alpha)))/nrow(D2_s1[D2_s1$R==0,]) +
        sum(ratio*(I(CFS2_sk<=theta)-m_sk_val[D2_sk$R==1]))/nrow(D2_sk[D2_sk$R==1,])
    }
    r_0 <- c(CFS1_sk, CFS2_sk)
    rhat_sk = try(uniroot(f_sk,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
    if(class(rhat_sk)!= "try-error"){
      rhat_sk=rhat_sk
    } else if (f_sk(max(r_0)/2)>0) {
      rhat_sk = try(uniroot(f_sk,c(0,max(r_0)/2),tol=1e-5)$root)
    } else {
      rhat_sk = try(uniroot(f_sk,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
    }
    rhat_ntgt[k-1] <- rhat_sk
    chi[k-1] <- abs(rhat_sk - rhat_tgSt)
    print("k:")
    print(k)
  }

  ### Calculating federated weights via optimization
  # [Optimized] Cache m_s1
  phi_s1 = function(theta) { # phi function of the target site
    m_s1_val <- m_s1(theta)
    ind1 <- which(D2$R==0 & D2$Site==1)
    ind2 <- which(D2$R==1 & D2$Site==1)
    phi <- rep(NA, length(D2$R))
    phi[ind1] <- m_s1_val[D2_s1$R==0]-(1-alpha)/mean(D2$Site==1 & D2$R==0)
    phi[ind2] <- eta_s1*(I(CFS2_s1<=theta)-m_s1_val[D2_s1$R==1])/mean(D2$Site==1 & D2$R==1)
    phi[-c(ind1,ind2)] <- 0
    return(matrix(phi, ncol=1))
  }
  # [Optimized] Direct glm in inner m_sk + cache all
  phi_k <- function(theta, k){
    D1_sk <- D1 %>% filter(Site==k)
    D2_sk <- D2 %>% filter(Site==k)
    D11_sk <- D11 %>% filter(Site==k)
    D21_sk <- D21 %>% filter(Site==k)

    m_sk <- function(theta) {
      y <- as.numeric(CFS1[[k-1]] <= theta)
      x <- D11_sk[, covar.names]
      fit <- glm(y ~ ., data = data.frame(y = y, x), family = binomial())
      p_CFS_sk <- as.numeric(predict(fit, newdata = D2_sk[, covar.names], type = "response"))
      return(p_CFS_sk)
    }

    m_s1_val <- m_s1(theta)
    m_sk_val <- m_sk(theta)
    ind1 <- which(D2$R==0 & D2$Site==1)
    ind2 <- which(D2$R==1 & D2$Site==k)
    phi <- rep(NA, length(D2$R))
    phi[ind1] <- (m_s1_val[D2_s1$R==0]-(1-alpha))/mean(D2$Site==1 & D2$R==0)
    phi[ind2] <- (omega[[k-1]]*(I(CFS2[[k-1]]<=theta)-m_sk_val[D2_sk$R==1]))/mean(D2$Site==k & D2$R==1)
    phi[-c(ind1,ind2)] <- 0

    inds1 <- which(D2$Site==1)
    phi[inds1] <- as.numeric(phi_s1(theta))[inds1]
    return(matrix(phi, ncol=1))
  }

  IF_s1 <- phi_s1(rhat_tgSt)
  IF_diff <- matrix(0, length(IF_s1), s-1)
  for(k in 2:s) {
    IF_diff[,(k-1)] <- phi_k(rhat_tgSt, k)
  }

  cvfit = cv.glmnet(x = IF_diff, y = IF_s1)
  lambda = cvfit$lambda.min
  fit = glmnet(x = IF_diff, y = IF_s1,
               penalty.factor = chi^2,
               intercept = FALSE,
               alpha = 1,
               lambda = lambda,
               lower.limits = 0,
               upper.limits = 1)
  w1 = coef(fit, s = lambda)[2:s]
  w2 = w1*(s-1)/s # to ensure at least 1/s weight assigned to target site

  fit3 = glmnet(x = IF_diff, y = IF_s1,
                penalty.factor = chi^2,
                intercept = FALSE,
                alpha = 1,
                lambda = lambda,
                lower.limits = 0,
                upper.limits = 1/s)
  w3 = coef(fit3, s = lambda)[2:s]

  ### Use federated weights
  rhat_fed1 = (1-sum(w1))*rhat_tgSt + sum(w1*rhat_ntgt)
  rhat_fed2 = (1-sum(w2))*rhat_tgSt + sum(w2*rhat_ntgt)
  rhat_fed3 = (1-sum(w3))*rhat_tgSt + sum(w3*rhat_ntgt)

  ### Use equal weights
  rhat_eqwt = mean(c(rhat_tgSt, rhat_ntgt))

  ###############################################################################
  ##############       Prediction on testing data point           ###############
  ###############################################################################
  ####data.test$R <- as.numeric(is.na(data.test[,Y.name])==F)
  ####data.test1 <- data.test %>% filter(R==1)
  ####data.test$R <- as.numeric(is.na(data.test[,Y.name])==F)
  ####data.test1 <- rbind(D11,D21)
  data.test1 <- rbind(D11_s1,D21_s1)####20240911

  if(conf.score=="ASR") {
    W=data.matrix(data.test1[,covar.names])
    #W.pred=data.matrix(data.test[data.test$R==0,covar.names])
    W.pred=data.matrix(Xtest2[,covar.names])
    Ytemp=c(data.test1[,Y.name])
    lambda_seq <- 10^seq(2, -2, by = -.1)
    mu <- glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    ridge_cv <- cv.glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    pred_y0 <- predict(mu, s=ridge_cv$lambda.min, W.pred)

    #PI.tgSt <- data.frame(lwr= -rhat_tgSt+pred_y0, upr=rhat_tgSt+pred_y0)
    PI.fed1 <- data.frame(lwr= -rhat_fed1+pred_y0, upr=rhat_fed1+pred_y0)
    PI.fed2 <- data.frame(lwr= -rhat_fed2+pred_y0, upr=rhat_fed2+pred_y0)
    PI.fed3 <- data.frame(lwr= -rhat_fed3+pred_y0, upr=rhat_fed3+pred_y0)
    PI.pool <- data.frame(lwr= -rhat_pool+pred_y0, upr=rhat_pool+pred_y0)
    PI.tgSt <- data.frame(lwr= -rhat_tgSt+pred_y0, upr=rhat_tgSt+pred_y0)
    PI.eqwt <- data.frame(lwr= -rhat_eqwt+pred_y0, upr=rhat_eqwt+pred_y0)
    #print(b,"finished")

  }
  if(conf.score=="localASR") {
    W=data.matrix(data.test1[,covar.names])
    #W.pred=data.matrix(data.test[data.test$R==0,covar.names])
    W.pred=data.matrix(Xtest2[,covar.names])
    Ytemp=data.test1[,Y.name]
    lambda_seq <- 10^seq(2, -2, by = -.1)
    mu <- glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    ridge_cv <- cv.glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    pred_y0 <- predict(mu, s=ridge_cv$lambda.min, W.pred)
    pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)

    ASR0 <- abs(pred_y1-Ytemp)
    Rho <- glmnet(W, ASR0, alpha=0, lambda=lambda_seq)
    Rho_cv <- cv.glmnet(W, ASR0, alpha=0, lambda=lambda_seq)
    rho0 <- predict(Rho, s=Rho_cv$lambda.min, W.pred)

    ####PI.tgSt <- data.frame(lwr= -rhat_tgSt*rho0+pred_y0, upr=rhat_tgSt*rho0+pred_y0)
     PI.fed1 <- data.frame(lwr= -rhat_fed1*rho0+pred_y0, upr=rhat_fed1*rho0+pred_y0)
    PI.fed2 <- data.frame(lwr= -rhat_fed2*rho0+pred_y0, upr=rhat_fed2*rho0+pred_y0)
    PI.fed3 <- data.frame(lwr= -rhat_fed3*rho0+pred_y0, upr=rhat_fed3*rho0+pred_y0)
    PI.pool <- data.frame(lwr= -rhat_pool*rho0+pred_y0, upr=rhat_pool*rho0+pred_y0)
    PI.tgSt <- data.frame(lwr= -rhat_tgSt*rho0+pred_y0, upr=rhat_tgSt*rho0+pred_y0)
    PI.eqwt <- data.frame(lwr= -rhat_eqwt*rho0+pred_y0, upr=rhat_eqwt*rho0+pred_y0)
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    fit.lw <- rq(rq.form, data=data.test1, tau=alpha/2)
    fit.up <- rq(rq.form, data=data.test1, tau=1-alpha/2)

    ####W.pred=data.test[data.test$R==0,]
    W.pred=Xtest2
    pred_y0.lw <- predict(fit.lw, W.pred)
    pred_y0.up <- predict(fit.up, W.pred)

    ####PI.tgSt <- data.frame(lwr= pred_y0.lw-rhat_tgSt, upr= pred_y0.up+rhat_tgSt)
    PI.fed1 <- data.frame(lwr= pred_y0.lw-rhat_fed1, upr= pred_y0.up+rhat_fed1)
    PI.fed2 <- data.frame(lwr= pred_y0.lw-rhat_fed2, upr= pred_y0.up+rhat_fed2)
    PI.fed3 <- data.frame(lwr= pred_y0.lw-rhat_fed3, upr= pred_y0.up+rhat_fed3)
    PI.pool <- data.frame(lwr= pred_y0.lw-rhat_pool, upr= pred_y0.up+rhat_pool)
    PI.tgSt <- data.frame(lwr= pred_y0.lw-rhat_tgSt, upr= pred_y0.up+rhat_tgSt)
    PI.eqwt <- data.frame(lwr= pred_y0.lw-rhat_eqwt, upr= pred_y0.up+rhat_eqwt)
  }
  ####colnames(PI.tgSt) <- c("lower", "upper")
  colnames(PI.fed1) <- colnames(PI.fed2) <- colnames(PI.fed3) <-
    colnames(PI.pool) <- colnames(PI.tgSt) <- colnames(PI.eqwt) <- c("lower", "upper")
    #PredInt[,1]=PI.tgSt$lower
    #PredInt[,2]=PI.tgSt$upper
    #YBleft <- rbind(YBleft,PredInt[,1])
    #YBright <- rbind(YBright,PredInt[,2])
    YBleft <- rbind(YBleft,PI.fed2$lower)
    YBright <- rbind(YBright,PI.fed2$upper)
  }

  #if bootstrap used, combine multiple CI to a single one
  #otherwise, a single CI is produced
  lefts <- NULL
  rights <- NULL
  Lens <- NULL
  if(B > 1){

    #ntest <- nrow(Xtest)####20250114
    ntest <- nrow(Xtest2)
    #ntest <- ncol(YBleft)####20250114
    for (j in 1:ntest) {

      Ys <- c(YBleft[,j],YBright[,j])
      Ws <- c(rep(1,nrow(YBleft)),rep(1,nrow(YBright)))
      Ss <- c(rep(1,nrow(YBleft)),rep(0,nrow(YBright)))
      #tol <- B*(1-gamma)
      tol <- B/2

      ind <- order(Ys)
      Ys <- Ys[ind]
      Ws <- Ws[ind]
      Ss <- Ss[ind]
      count <- 0
      leftpoints <- NULL
      rightpoints <- NULL

      for (i in 1:length(Ys)) {
        if(Ss[i]==1){
          count <- count+Ws[i]
          if((count>tol)&(count-Ws[i]<=tol)){
            leftpoints <- c(leftpoints,Ys[i])
          }
        }else if(Ss[i]==0){
          if((count>tol)&(count-Ws[i]<=tol)){
            rightpoints <- c(rightpoints,Ys[i])
          }
          count <- count-Ws[i]
        }
      }
      lefts[[j]] <- leftpoints
      rights[[j]] <- rightpoints
      Lens[j] <- sum(rightpoints-leftpoints)
    }
  } else if (B == 1) {                 # B==1 FIX: single split, works for any method
    ntest <- ncol(YBleft)
    for (j in 1:ntest) {
      lefts[[j]]  <- YBleft[1, j]
      rights[[j]] <- YBright[1, j]
      Lens[j]     <- YBright[1, j] - YBleft[1, j]
    }
  }
  ####restore later
  return(list(lefts=lefts, rights=rights, Lens=Lens, Nts=Nts,r0=r0))

}

ASR <- function(W, Ytemp, W.pred, Y.pred) {
  lambda_seq <- 10^seq(2, -2, by = -.1)
  mu <- glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
  ridge_cv <- cv.glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
  pred_y <- predict(mu, s=ridge_cv$lambda.min, W.pred)
  pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)

  CFS1 <- as.numeric(abs(pred_y1-Ytemp))
  CFS2 <- as.numeric(abs(pred_y-Y.pred))
  return(CFS=list(CFS1=CFS1, CFS2=CFS2))
}

localASR <- function(W, Ytemp, W.pred, Y.pred) {
  lambda_seq <- 10^seq(2, -2, by = -.1)
  mu <- glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
  ridge_cv <- cv.glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)

  pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)
  pred_y <- predict(mu, s=ridge_cv$lambda.min, W.pred)

  ASR1 <- abs(pred_y1-Ytemp)
  ASR <- abs(pred_y-Y.pred)

  Rho <- glmnet(W, ASR1, alpha=0, lambda=lambda_seq)
  Rho_cv <- cv.glmnet(W, ASR1, alpha=0, lambda=lambda_seq)
  rho1 <- predict(Rho, s=Rho_cv$lambda.min, W)
  rho <- predict(Rho, s=Rho_cv$lambda.min, W.pred)

  CFS1 <- as.numeric(ASR1/rho1)
  CFS2 <- as.numeric(ASR/rho)
  return(CFS=list(CFS1=CFS1, CFS2=CFS2))
}

CQR <- function(data, Ytemp, data.pred, Y.pred, rq.form, alpha=0.1) {
  fit.lw <- rq(rq.form, data=data, tau=alpha/2)
  fit.up <- rq(rq.form, data=data, tau=1-alpha/2)

  pred_y1.lw <- predict(fit.lw, data)
  pred_y1.up <- predict(fit.up, data)

  pred_y.lw <- predict(fit.lw, data.pred)
  pred_y.up <- predict(fit.up, data.pred)

  CFS1 <- as.numeric(apply(cbind(pred_y1.lw-Ytemp, Ytemp-pred_y1.up), 1, max))
  CFS2 <- as.numeric(apply(cbind(pred_y.lw-Y.pred, Y.pred-pred_y.up), 1, max))
  return(CFS=list(CFS1=CFS1, CFS2=CFS2))
}

MuSCI_ConformalCI <- function(X,Y,T,Deta,Xtest,method,B,trainprop,psfun,psparams,outfun,outparams,quantiles,type,side,weighted,
####                  conf.score="ASR",
                  conf.score,
                  ps.Library="SL.glm",
                  m.Library="SL.glm",
                  p1,#0.5
                  seed=10012) {
  YBleft <- NULL
  YBright <- NULL
  Nts <- NULL
  r0<-NULL
  r1<-NULL
#random sample from binary varibale(A) B times
  for (b in 1:B) {
  n <- length(Y)
  #split samples for training and calibration
  YY <- Y
  YY[!T] <- NA
  T <- as.numeric(!is.na(YY))
  res <- SampleSplit(YY,trainprop)
  #res$trainid: index for all training
  #res$trainid1: index for all training with T=1
  #!!!need to select T=1 first and then select trainid1
  #res$trainid0: index for all training with T=0

  #estimate the propensity score
  Xtrain <- X[res$trainid, ,drop=FALSE]
  Ttrain <- T[res$trainid]
  PSFit_T <- PropScore(Ttrain,Xtrain,psfun, psparams)

  #estimate weights for testing data
  ps_T_test <- PSFit_T(Xtest)
  Xtest <- as.matrix(Xtest)
  ps_D_test <- exp(cbind(1,Xtest)%*%Deta)/(1+exp(cbind(1,Xtest)%*%Deta))
  ps_test <- ps_T_test*(1-ps_T_test)/(ps_T_test*(1-ps_D_test)+ps_D_test*(1-ps_T_test))


  #estimate the propensity score for binary A
  tp_T <- PSFit_T(X)
  tp_D <- exp(cbind(1,X)%*%Deta)/(1+exp(cbind(1,X)%*%Deta))
  tp_A <- tp_D*(1-tp_T)/(tp_D*(1-tp_T)+tp_T*(1-tp_D))


    A <- rbinom(n = n, size = 1, prob=tp_A)
    ##print("Len:")
    ##print(length(T)[1])
    ##print(n)
    ##print(length(A))
    O <- I(T==A)

    Nts[b] <- sum(O)
    #quantile regression fitting after selecting A=T
    X1 <- X[res$trainid, ,drop=FALSE]
    Y1 <- Y[res$trainid]
    #select O=1 from training data
    ind <- which(O[res$trainid] == 1)
    Xtrain <- X1[ind, ,drop=FALSE]
    Ytrain <- Y1[ind]
    ####OutputFit <- QuantileReg(Ytrain,Xtrain,outfun,outparams,quantiles)
    ####Yhat_test <- OutputFit(Xtest)
    ###sample and fit for train

    ####sample for val
    if(method == "Sweight" | method == "BootSweight"){
      #data for validation
      X0 <- X[-res$trainid, ,drop=FALSE]
      Y0 <- Y[-res$trainid]
      ind <- which(O[-res$trainid] == 1)
      Xval <- X0[ind, ,drop=FALSE]
      Yval <- Y0[ind]
      ####ps_T_val <- PSFit_T(as.data.frame(Xval))
      ####Yhat_val <- OutputFit(Xval)

      Xval <- as.matrix(Xval)
      ####ps_D_val <- exp(cbind(1,Xval)%*%Deta)/(1+exp(cbind(1,Xval)%*%Deta))
      ####ps_val <- ps_T_val*(1-ps_T_val)/(ps_T_val*(1-ps_D_val)+ps_D_val*(1-ps_T_val))

      ####subsample wttl, wttl is for both train and val
      ####wttl0 = wttl[-res$trainid]
      ####wttl_val = wttl0[ind]
    }
    #else


  ###############################################################################
  ##################            Data pre-processing           ###################
  ###############################################################################
  ####colnames(data.train)[colnames(data.train)==T.name] <- "Site"
  ####site.names <- unique(data.train$Site) # extract names of sites
  ####s <- length(site.names) # number of sites
  ####s=1
  ####non.tgt.site.names <- site.names[site.names%in%tgt.site.name==F] # names of non-target sites

  ####data.train$Site[data.train$Site==tgt.site.name] <- 999
  ####Xtrain$Site <- 999
  ###for(k in 2:s) {
  ###  data.train$Site[data.train$Site==non.tgt.site.names[k-1]] <- k*999
  ###}
  ####Xtrain$Site <- Xtrain$Site/999
  ####covar.names <- names(Xtrain)[names(Xtrain)%in%c(Y.name, "Site")==F]
  D1 <- D2 <- data.frame()
  set.seed(seed)

  D11=data.frame(Xtrain)
  D21=data.frame(Xval)
  Xtest2=data.frame(Xtest)
  colnames(D11)=colnames(D21)=colnames(Xtest2)=paste0("cn",1:ncol(D11))
  ###############################################################################
  ##################            Pooled sample method           ##################
  ###############################################################################
  ####D1$R <- as.numeric(is.na(D1[,Y.name])==F) # define the missingness indicator
  ####D2$R <- as.numeric(is.na(D2[,Y.name])==F)
  ####D11 <- D1[D1$R==1,]
  ####D21 <- D2[D2$R==1,]

  ####y = D1$R
  ####x = D1[,covar.names]

  y=c(rep(1,nrow(rbind(D11,D21))),rep(0,nrow(Xtest2)))
  x=rbind(D11,D21,Xtest2)
  ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
  pred = as.numeric(predict(ps.fit, D21)$library.predict)
  prob = pmax(pmin(pred, 0.9), 0.1)
  eta = (1-prob)/prob

  D11all=cbind(D11,Ytrain)
    D21all=cbind(D21,Yval)
    colnames(D11all)[ncol(D11all)]="Yname"
    colnames(D21all)[ncol(D21all)]="Yname"
    Y.name="Yname"
    D11=D11all
    D21=D21all
    res2 <- SampleSplit(Xtest2[,1],trainprop)
    D10=data.frame(cbind(Xtest2[res2$trainid,],rep(NA,nrow(Xtest[res2$trainid,]))))
    colnames(D10)=colnames(D11)
    D1=rbind(D11,D10)
    D20=data.frame(cbind(Xtest2[-res2$trainid,],rep(NA,nrow(Xtest[-res2$trainid,]))))
    colnames(D20)=colnames(D21)
    D2=rbind(D21,D20)
    colnames(D1)=colnames(D11)
    colnames(D2)=colnames(D21)
    D1$R=D2$R=0
    D1$R[1:nrow(D11)]=1
    D2$R[1:nrow(D21)]=1
    covar.names=colnames(D21)[1:(ncol(D21)-1)]

  if(conf.score=="ASR") {
    ##print("Y")
    ##print(head(D11[,Y.name]))
    CFS <- ASR(W=data.matrix(D11[,covar.names]),
               Ytemp=D11[,Y.name],
               W.pred=data.matrix(D21[,covar.names]),
               Y.pred=D21[,Y.name])
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  ##print("flag1")
  if(conf.score=="localASR") {
    CFS <- localASR(W=data.matrix(D11[,covar.names]),
                    Ytemp=D11[,Y.name],
                    W.pred=data.matrix(D21[,covar.names]),
                    Y.pred=D21[,Y.name])
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    CFS <- CQR(data=D11, Ytemp=D11[,Y.name],
               data.pred=D21, Y.pred=D21[,Y.name],
               rq.form=rq.form, alpha=alpha)
    CFS1 <- CFS[[1]]; CFS2 <- CFS[[2]]
  }

  # [Optimized] Direct glm + cached m(theta)
  m <- function(theta) {
    y <- as.numeric(CFS1 <= theta)
    x <- D11[, covar.names]
    fit <- glm(y ~ ., data = data.frame(y = y, x), family = binomial())
    p_CFS <- as.numeric(predict(fit, newdata = D2[, covar.names], type = "response"))
    return(p_CFS)
  }
  f <- function(theta){
    m_val <- m(theta)
    sum(m_val[D2$R==0]-(1-alpha))/nrow(D2[D2$R==0,]) +
       sum(eta*(I(CFS2<=theta)-m_val[D2$R==1]))/nrow(D2[D2$R==1,])
  }
  r_0 <- c(CFS1, CFS2)
  rhat_pool = try(uniroot(f,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
  if(class(rhat_pool)!= "try-error"){
    rhat_pool = rhat_pool
  } else if (f(max(r_0)/2)>0) {
    rhat_pool = try(uniroot(f,c(0,max(r_0)/2),tol=1e-5)$root)
  } else {
    rhat_pool = try(uniroot(f,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
  }
  ###############################################################################
  ##################         Use target site data only         ##################
  ###############################################################################
  ####D1_s1 <- D1[D1$Site==1, ]
  ####D2_s1 <- D2[D2$Site==1, ]
  ####D11_s1 <- D1_s1[D1_s1$R==1, ]
  ####D21_s1 <- D2_s1[D2_s1$R==1, ]
  D1_s1 <- D1
  D2_s1 <- D2
  D11_s1 <- D11
  D21_s1 <- D21

  y = D1_s1$R
  x = D1_s1[,covar.names]
  ps.fit = SuperLearner(Y=y, X=x, family=binomial(), SL.library=ps.Library)
  pred = as.numeric(predict(ps.fit, D21_s1[,covar.names])$library.predict)
  prob = pmax(pmin(pred, 0.9), 0.1)
  eta_s1 = (1-prob)/prob

  if(conf.score=="ASR") {
    CFS <- ASR(W=data.matrix(D11_s1[,covar.names]),
               Ytemp=D11_s1[,Y.name],
               W.pred=data.matrix(D21_s1[,covar.names]),
               Y.pred=D21_s1[,Y.name])
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  if(conf.score=="localASR") {
    CFS <- localASR(W=data.matrix(D11_s1[,covar.names]),
                    Ytemp=D11_s1[,Y.name],
                    W.pred=data.matrix(D21_s1[,covar.names]),
                    Y.pred=D21_s1[,Y.name])
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    CFS <- CQR(data=D11_s1, Ytemp=D11_s1[,Y.name],
               data.pred=D21_s1, Y.pred=D21_s1[,Y.name],
               rq.form=rq.form, alpha=alpha)
    CFS1_s1 <- CFS[[1]]; CFS2_s1 <- CFS[[2]]
  }

  # [Optimized] Direct glm + cached m_s1
  m_s1 <- function(theta) {
    y <- as.numeric(CFS1_s1 <= theta)
    x <- D11_s1[, covar.names]
    fit <- glm(y ~ ., data = data.frame(y = y, x), family = binomial())
    p_CFS <- as.numeric(predict(fit, newdata = D2_s1[, covar.names], type = "response"))
    return(p_CFS)
  }

  f_s1 <- function(theta){
    m_s1_val <- m_s1(theta)
    sum(m_s1_val[D2_s1$R==0]-(1-alpha))/nrow(D2_s1[D2_s1$R==0,])  +
       sum(eta_s1*(I(CFS2_s1<=theta)-m_s1_val[D2_s1$R==1]))/nrow(D2_s1[D2_s1$R==1,])
  }
  r_0 <- c(CFS1_s1, CFS2_s1)
  rhat_tgSt = try(uniroot(f_s1,c(1,max(r_0)/2),extendInt="yes",tol=1e-5)$root)
  if(class(rhat_tgSt)!= "try-error"){
    rhat_tgSt = rhat_tgSt
  } else if (f_s1(max(r_0)/2)>0) {
    rhat_tgSt = try(uniroot(f_s1,c(0,max(r_0)/2),tol=1e-5)$root)
  } else {
    rhat_tgSt = try(uniroot(f_s1,c(max(r_0)/2,max(r_0)+0.5),tol=1e-5)$root)
  }
  ####
  r0=c(r0,rhat_tgSt)
print("flag")
  ###############################################################################
  ##############       Prediction on testing data point           ###############
  ###############################################################################
  ####data.test$R <- as.numeric(is.na(data.test[,Y.name])==F)
  ####data.test1 <- data.test %>% filter(R==1)
  ####data.test$R <- as.numeric(is.na(data.test[,Y.name])==F)
  data.test1 <- rbind(D11,D21)

  if(conf.score=="ASR") {
    W=data.matrix(data.test1[,covar.names])
    #W.pred=data.matrix(data.test[data.test$R==0,covar.names])
    W.pred=data.matrix(Xtest2[,covar.names])
    Ytemp=c(data.test1[,Y.name])
    lambda_seq <- 10^seq(2, -2, by = -.1)
    mu <- glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    ridge_cv <- cv.glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    pred_y0 <- predict(mu, s=ridge_cv$lambda.min, W.pred)

    PI.tgSt <- data.frame(lwr= -rhat_tgSt+pred_y0, upr=rhat_tgSt+pred_y0)
    #print(b,"finished")

  }
  if(conf.score=="localASR") {
    W=data.matrix(data.test1[,covar.names])
    #W.pred=data.matrix(data.test[data.test$R==0,covar.names])
    W.pred=data.matrix(Xtest2[,covar.names])
    Ytemp=data.test1[,Y.name]
    lambda_seq <- 10^seq(2, -2, by = -.1)
    mu <- glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    ridge_cv <- cv.glmnet(W, Ytemp, alpha=0, lambda=lambda_seq)
    pred_y0 <- predict(mu, s=ridge_cv$lambda.min, W.pred)
    pred_y1 <- predict(mu, s=ridge_cv$lambda.min, W)

    ASR0 <- abs(pred_y1-Ytemp)
    Rho <- glmnet(W, ASR0, alpha=0, lambda=lambda_seq)
    Rho_cv <- cv.glmnet(W, ASR0, alpha=0, lambda=lambda_seq)
    rho0 <- predict(Rho, s=Rho_cv$lambda.min, W.pred)

    PI.tgSt <- data.frame(lwr= -rhat_tgSt*rho0+pred_y0, upr=rhat_tgSt*rho0+pred_y0)

  }
  if(conf.score=="CQR") {
    rq.form <- as.formula(paste(Y.name, "~", paste(covar.names, collapse="+")))
    fit.lw <- rq(rq.form, data=data.test1, tau=alpha/2)
    fit.up <- rq(rq.form, data=data.test1, tau=1-alpha/2)

    ####W.pred=data.test[data.test$R==0,]
    W.pred=Xtest2
    pred_y0.lw <- predict(fit.lw, W.pred)
    pred_y0.up <- predict(fit.up, W.pred)

    PI.tgSt <- data.frame(lwr= pred_y0.lw-rhat_tgSt, upr= pred_y0.up+rhat_tgSt)

  }
  colnames(PI.tgSt) <- c("lower", "upper")
    #PredInt[,1]=PI.tgSt$lower
    #PredInt[,2]=PI.tgSt$upper
    #YBleft <- rbind(YBleft,PredInt[,1])
    #YBright <- rbind(YBright,PredInt[,2])
    YBleft <- rbind(YBleft,PI.tgSt$lower)
    YBright <- rbind(YBright,PI.tgSt$upper)
  }

  #if bootstrap used, combine multiple CI to a single one
  #otherwise, a single CI is produced
  lefts <- NULL
  rights <- NULL
  Lens <- NULL
  if(B > 1){

    ntest <- nrow(Xtest)
    for (j in 1:ntest) {

      Ys <- c(YBleft[,j],YBright[,j])
      Ws <- c(rep(1,nrow(YBleft)),rep(1,nrow(YBright)))
      Ss <- c(rep(1,nrow(YBleft)),rep(0,nrow(YBright)))
      #tol <- B*(1-gamma)
      tol <- B/2

      ind <- order(Ys)
      Ys <- Ys[ind]
      Ws <- Ws[ind]
      Ss <- Ss[ind]
      count <- 0
      leftpoints <- NULL
      rightpoints <- NULL

      for (i in 1:length(Ys)) {
        if(Ss[i]==1){
          count <- count+Ws[i]
          if((count>tol)&(count-Ws[i]<=tol)){
            leftpoints <- c(leftpoints,Ys[i])
          }
        }else if(Ss[i]==0){
          if((count>tol)&(count-Ws[i]<=tol)){
            rightpoints <- c(rightpoints,Ys[i])
          }
          count <- count-Ws[i]
        }
      }
      lefts[[j]] <- leftpoints
      rights[[j]] <- rightpoints
      Lens[j] <- sum(rightpoints-leftpoints)
    }
  } else if (B == 1) {                 # B==1 FIX: single split, works for any method
    ntest <- ncol(YBleft)
    for (j in 1:ntest) {
      lefts[[j]]  <- YBleft[1, j]
      rights[[j]] <- YBright[1, j]
      Lens[j]     <- YBright[1, j] - YBleft[1, j]
    }
  }
  ####restore later
  return(list(lefts=lefts, rights=rights, Lens=Lens, Nts=Nts,r0=r0))

}

Estimate_omega_np = function(x, x_target, x.pred, x_target.pred) {
  x.all = rbind(x, x_target)
  z.all = c(rep(1,nrow(x)),rep(0,nrow(x_target)))

  fit <- cv.glmnet(x=x.all, y=z.all, nfolds=5, family='binomial')
  src.predict = predict(fit, x.pred, type="response", lambda=fit$lambda.min)

  omega = (1-src.predict)/src.predict # *nrow(x.pred)/nrow(x_target.pred)
  omega = pmax(pmin(omega, 20), 0.05)
  return(omega)
}

genY <- function(X){
  n <- nrow(X)
  2 / (1 + exp(-12 * (X[, 1] - 0.5))) * 2 / (1 + exp(-12 * (X[, 2] - 0.5))) + rnorm(n)
}

genY_singlestage <- function(X,T){
  n <- length(T)
  Y <- 1+X[,1]-X[,2]+X[,3]^{3}+exp(X[,4])+T*(3-5*X[,1]+2*X[,2]-3*X[,3]+X[,4])+(1+T)*(1+X[,1]+X[,2]+X[,3]+X[,4])*rnorm(n)
  return(Y)
}

genY_twostage <- function(X1,X2,T1,T2){
  n <- length(T1)
  Y <- 1+X1+T1*(1-3*(X1-0.2)^{2})+X2+T2*(1-5*(X2-0.4)^{2})+(1+0.5*T1-T1*X1+0.5*T2-T2*X2)*rnorm(n,0,sd=0.4)
  return(Y)
  }

ConformalCI <- function(X,Y,T,Deta,Xtest,method,B=1,trainprop,psfun,psparams,outfun,outparams,quantiles,type,side,weighted,alpha,wthigh,wtlow,wttl,wttl_test){
####added wttl as the weight for covariate shift

  YBleft <- NULL
  YBright <- NULL
  Nts <- NULL
  #random sample from binary varibale(A) B times
  for (b in 1:B) {
  n <- length(Y)
  #split samples for training and calibration
  YY <- Y
  YY[!T] <- NA
  T <- as.numeric(!is.na(YY))
  res <- SampleSplit(YY,trainprop)
  #res$trainid: index for all training
  #res$trainid1: index for all training with T=1
  #!!!need to select T=1 first and then select trainid1
  #res$trainid0: index for all training with T=0

  #estimate the propensity score
  Xtrain <- X[res$trainid, ,drop=FALSE]
  Ttrain <- T[res$trainid]
  PSFit_T <- PropScore(Ttrain,Xtrain,psfun, psparams)

  #estimate weights for testing data
  ps_T_test <- PSFit_T(Xtest)
  Xtest <- as.matrix(Xtest)
  ps_D_test <- exp(cbind(1,Xtest)%*%Deta)/(1+exp(cbind(1,Xtest)%*%Deta))
  ps_test <- ps_T_test*(1-ps_T_test)/(ps_T_test*(1-ps_D_test)+ps_D_test*(1-ps_T_test))


  #estimate the propensity score for binary A
  tp_T <- PSFit_T(X)
  tp_D <- exp(cbind(1,X)%*%Deta)/(1+exp(cbind(1,X)%*%Deta))
  tp_A <- tp_D*(1-tp_T)/(tp_D*(1-tp_T)+tp_T*(1-tp_D))


    A <- rbinom(n = n, size = 1, prob=tp_A)
    O <- I(T==A)
    Nts[b] <- sum(O)
    #quantile regression fitting after selecting A=T
    X1 <- X[res$trainid, ,drop=FALSE]
    Y1 <- Y[res$trainid]
    #select O=1 from training data
    ind <- which(O[res$trainid] == 1)
    Xtrain <- X1[ind, ,drop=FALSE]
    Ytrain <- Y1[ind]
    OutputFit <- QuantileReg(Ytrain,Xtrain,outfun,outparams,quantiles)
    Yhat_test <- OutputFit(Xtest)

    if(method == "Sweight" | method == "BootSweight"){
      #data for validation
      X0 <- X[-res$trainid, ,drop=FALSE]
      Y0 <- Y[-res$trainid]
      ind <- which(O[-res$trainid] == 1)
      Xval <- X0[ind, ,drop=FALSE]
      Yval <- Y0[ind]
      ps_T_val <- PSFit_T(as.data.frame(Xval))
      Yhat_val <- OutputFit(Xval)

      Xval <- as.matrix(Xval)
      ps_D_val <- exp(cbind(1,Xval)%*%Deta)/(1+exp(cbind(1,Xval)%*%Deta))
      ps_val <- ps_T_val*(1-ps_T_val)/(ps_T_val*(1-ps_D_val)+ps_D_val*(1-ps_T_val))

      ####subsample wttl, wttl is for both train and val
      wttl0 = wttl[-res$trainid]
      wttl_val = wttl0[ind]


      }else if(method == "Dweight" | method== "BootDweight"){
      #regression for validation
      Xval <- X[-res$trainid, ,drop=FALSE]
      Yval <- Y[-res$trainid]
      ps_T_val <- PSFit_T(as.data.frame(Xval))
      Yhat_val <- OutputFit(Xval)

      Xval <- as.matrix(Xval)
      ps_D_val <- exp(cbind(1,Xval)%*%Deta)/(1+exp(cbind(1,Xval)%*%Deta))
      ps_val <- ps_T_val*(1-ps_T_val)/(ps_T_val*(1-ps_D_val)+ps_D_val*(1-ps_T_val))
      ps_A_val <- ps_D_val*(1-ps_T_val)/(ps_D_val*(1-ps_T_val)+ps_T_val*(1-ps_D_val))
      ps_A_val[which(T[-res$trainid]==0)] <- 1-ps_A_val[which(T[-res$trainid]==0)]
      ps_val <- ps_val/ps_A_val

      #### wttl is for both train and val
      wttl_val= wttl[-res$trainid]
    }

    Yscore <- conformalScore(Yval, Yhat_val, type, side)
    #construct the confidence interval by weighted CP
    if(length(weighted)>1){
      wt <- weighted[1:n]
      if(method == "Sweight" | method == "BootSweight"){
        wt <- wt[-res$trainid]
        ind <- which(O[-res$trainid] == 1)
        wt <- wt[ind]
      }else if(method == "Dweight" | method == "BootDweight"){
        wt <- wt[-res$trainid]
        }
      wt_test <- weighted[-(1:n)]
      PredInt<- predict.conformalSplit(Yhat_test, Yscore,wt,wt_test,alpha,wthigh,wtlow,type,side)
    }else if(weighted==TRUE){
      ####wt <- 1/ps_val
      ####add wttl, wttl_test
      wt <- (1/ps_val)*wttl_val
      wt_test <- (1/ps_test)*wttl_test
      PredInt<- predict.conformalSplit(Yhat_test, Yscore,wt,wt_test,alpha,wthigh,wtlow,type,side)
    }else if(weighted==FALSE){
      wt <- rep(1,length(ps_val))
      wt_test <- rep(1,length(ps_test))
      PredInt<- predict.conformalSplit(Yhat_test, Yscore,wt,wt_test,alpha,wthigh,wtlow,type,side)
    }

    YBleft <- rbind(YBleft,PredInt[,1])
    YBright <- rbind(YBright,PredInt[,2])
  }

  #if bootstrap used, combine multiple CI to a single one
  #otherwise, a single CI is produced
  lefts <- NULL
  rights <- NULL
  Lens <- NULL
  if(B > 1){

    ntest <- nrow(Xtest)
    for (j in 1:ntest) {

      Ys <- c(YBleft[,j],YBright[,j])
      Ws <- c(rep(1,nrow(YBleft)),rep(1,nrow(YBright)))
      Ss <- c(rep(1,nrow(YBleft)),rep(0,nrow(YBright)))
      #tol <- B*(1-gamma)
      tol <- B/2

      ind <- order(Ys)
      Ys <- Ys[ind]
      Ws <- Ws[ind]
      Ss <- Ss[ind]
      count <- 0
      leftpoints <- NULL
      rightpoints <- NULL

      for (i in 1:length(Ys)) {
        if(Ss[i]==1){
          count <- count+Ws[i]
          if((count>tol)&(count-Ws[i]<=tol)){
            leftpoints <- c(leftpoints,Ys[i])
          }
        }else if(Ss[i]==0){
          if((count>tol)&(count-Ws[i]<=tol)){
            rightpoints <- c(rightpoints,Ys[i])
          }
          count <- count-Ws[i]
        }
      }
      lefts[[j]] <- leftpoints
      rights[[j]] <- rightpoints
      Lens[j] <- sum(rightpoints-leftpoints)
    }
  } else if (B == 1) {                 # B==1 FIX: single split, works for any method
    ntest <- ncol(YBleft)
    for (j in 1:ntest) {
      lefts[[j]]  <- YBleft[1, j]
      rights[[j]] <- YBright[1, j]
      Lens[j]     <- YBright[1, j] - YBleft[1, j]
    }
  }

  return(list(lefts=lefts, rights=rights, Lens=Lens, Nts=Nts))
}

delta.tilde.prob <- function(x, alphas) {
  # x: covariate vecter or matrix including intercep 1
  # return: the probability of \widetilde{\delta}=1 given x
  x.alpha <- x %*% alphas
  prob <- exp(x.alpha)/(1 + exp(x.alpha))
  return(prob)
}

cost2 <- function(alphas, x.in.rct, x.in.rwe, N) {
  # the cost function of the first parametric weighting method
  part1 <- sum(x.in.rct %*% as.matrix(alphas, length(alphas), 1))
  part2 <- mean(log(1 + exp(x.in.rwe %*% as.matrix(alphas, length(alphas), 1))))*N
  return(part2-part1)
}

cost <- function(alphas, x.in.rct, left.hand, N) {
  # the cost function of the second parametric weighting method
  right.hand <- apply(x.in.rct / c(delta.tilde.prob(x.in.rct, alphas)), 2, sum)/N
  diff <- right.hand - left.hand
  return(sqrt(mean(diff**2)))
}

entr_balance <- function(tol, x.in.rct, x.in.rwe) {
  # the optimization problem for the nonparametric weighting method
  n <- dim(x.in.rct)[1]
  w <- Variable(n)
  p <- dim(x.in.rct)[2]-1
  # entr: is the element-wise entropy atom
  objective <- Maximize(sum(entr(w)))
  constraints <- list(w >= 0, sum(w) == 1)
  for (k in 1:p) {
    constraints <- c(constraints, abs(sum(w*x.in.rct[,k+1]) - mean(x.in.rwe[,k+1])) <= tol*sd(x.in.rwe[,k+1]))
    constraints <- c(constraints, abs(sum(w*(x.in.rct[,k+1]**2)) - mean(x.in.rwe[,k+1]**2)) <= tol*sd(x.in.rwe[,k+1]**2)) ## added for the misspecify
  }
  prob <- Problem(objective, constraints)
  w.tilde <- rep(0, n)
  tryCatch(
    expr = {
      result <- suppressWarnings(solve(prob))
      w.tilde <- result$getValue(w)
      if (is.na(w.tilde)[1]) {
        w.tilde <- rep(0,n)
      }
    },
    error = function(e) {
      # print("error")
    }
  )
  return(w.tilde)
}

weight_tuned_entr <- function(tols, x.in.rct, x.in.rwe, prop=0.5, smps=10) {
  # Evaluate covariate balance with bootstrap samples to choose tuning parameters:
  #   the degree of approximate balance sigma_k
  n <- dim(x.in.rct)[1]
  sds <- apply(x.in.rwe, 2, sd)
  means <- apply(x.in.rwe, 2, mean)
  sds2 <- apply(x.in.rwe**2, 2, sd)
  moment2 <- apply(x.in.rwe**2, 2, mean)
  cov.diff.bars <- c()
  w.tildes <- matrix(NA, n, length(tols))
  for (i in 1:length(tols)) {
    tol <- tols[i]
    w.tilde <- entr_balance(tol, x.in.rct, x.in.rwe)
    w.tilde[w.tilde < 0] <- 0 ## the optimization results sometimes return negative weights, so clip them to 0
    if (sum(w.tilde)==0) {
      cov.diff.bar <- 1e+10
    } else {
      w.tildes[,i] <- w.tilde
      cov.diffs <- c()
      for (s in 1:smps) {
        boot.ind <- sample.int(n, round(prop*n), replace=TRUE)
        x.in.rct.boot <- x.in.rct[boot.ind,]
        w.tilde.boot <- w.tilde[boot.ind]
        cov.diff <- (apply(x.in.rct.boot*w.tilde.boot, 2, sum)/sum(w.tilde.boot) -  means)[-1] / sds[-1] ## remember check this again, should also add 2nd moment diff?
        cov.diff.moment2 <- (apply(x.in.rct.boot**2*w.tilde.boot, 2, sum)/sum(w.tilde.boot) -  moment2)[-1] / sds2[-1]
        cov.diffs <- c(cov.diffs, sum((c(cov.diff, cov.diff.moment2))**2)) ## L2 measure
      }
      cov.diff.bar <- mean(cov.diffs)
    }
    cov.diff.bars <- c(cov.diff.bars, cov.diff.bar)
  }
  tuned.tol <- tols[which.min(cov.diff.bars)]
  if (min(cov.diff.bars) < 1e+10) {
    return(w.tildes[, which.min(cov.diff.bars)])
  } else {
    stop('************ optimization not feasible *****************')
  }
}

learn_weights <- function(x.in.rct, x.in.rwe, w.method,  N) {
  # w.method: 1, 2, 3, 4
  n <- dim(x.in.rct)[1]
  p <- dim(x.in.rct)[2] - 1
  if (w.method == 1) {
    # unweighted
    w.tilde <- rep(1, n)
    w.tilde.rwe <- rep(1, dim(x.in.rwe)[1])
  } else if (w.method == 2) {
    # the first parametric weighting method

      opt2 <- optim(rep(0,p+1), cost2, x.in.rct=x.in.rct, x.in.rwe=x.in.rwe, N=N, method = "Nelder-Mead", control = list(maxit=10000))
      (alphas.hat2 <- opt2$par)
      w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))
      w.tilde.rwe <- 1 + 1/(exp(x.in.rwe %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))

  } else if (w.method == 3) {
    # the second parametric weighting method
    left.hand <- apply(x.in.rwe, 2, mean)

      opt <- optim(rep(0, p+1), cost, x.in.rct=x.in.rct, left.hand=left.hand, N=N, method = "Nelder-Mead", control = list(maxit=10000))
      (alphas.hat3 <- opt$par)
      w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat3, length(alphas.hat3), 1)))
      w.tilde.rwe <- 1 + 1/(exp(x.in.rwe %*% as.matrix(alphas.hat3, length(alphas.hat3), 1)))

  } else if (w.method == 4) {
    # nonparametric weighting
    w.tilde <- weight_tuned_entr(tols=c(0.000, 0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500, 1.000),
                                 x.in.rct, x.in.rwe)
  } else if (w.method == 5) {
    # the first parametric weighting method with true alphas

      #opt2 <- optim(rep(0,p+1), cost2, x.in.rct=x.in.rct, x.in.rwe=x.in.rwe, N=N, method = "Nelder-Mead", control = list(maxit=10000))
      #(alphas.hat2 <- opt2$par)
    (alphas.hat2 <- c(alphas.true))
      w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))
      w.tilde.rwe <- 1 + 1/(exp(x.in.rwe %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))

  } else {
    stop("Wrong input for weighting method")
  }
  ####w.tilde
  return(list(w.tilde.rct=w.tilde, w.tilde.rwe=w.tilde.rwe))
}

find_inds <- function(a, b){
    n <- length(a)
    b <- b - 1e-12
    ## n + 1 - rank(-c(a, b), ties.method = "first")[-(1:n)] + rank(-b, ties.method = "first")
    rank(c(a, b), ties.method = "first")[-(1:n)] - rank(b, ties.method = "first") + 1
}

weightedConformalCutoff <- function(score, weight, qt){
    ord <- order(score)
    weight <- weight[ord]
    score <- score[ord]
    cw <- cumsum(weight)
    inds <- find_inds(cw, qt)
    cutoff <- score[inds]
    return(cutoff)
}

conformalScore <- function(Y, Yhat, type, side){
    if (is.vector(Y) || (is.matrix(Y) && ncol(Y) == 1)){
        if (type == "CQR" && side == "two"){
            score <- pmax(Yhat[, 1] - Y, Y - Yhat[, 2])
        } else if (type == "mean" && side == "two"){
            score <- abs(Yhat - Y)
        } else if (side == "above"){
            score <- Y - Yhat
        } else if (side == "below"){
            score <- Yhat - Y
        }
    } else if (is.matrix(Y) && ncol(Y) == 2){
        score <- pmax(Yhat[, 1] - Y[, 1], Y[, 2] - Yhat[, 2])
    }
    return(score)
}

conformalScoreClass <- function(Y, phat, type, wt){
    ncl <- nlevels(Y)
    if (type == "weighted"){
        phat <- row_quo(phat, wt)
    }
    Yid <- as.numeric(Y)
    score <- phat[cbind(1:length(Y), Yid)]
    return(score)
}

censoring <- function(x, high = 20, low = 0.05){
    pmin(pmax(x, low), high)
}

predict.conformalSplit <- function(Yhat_test,Yscore, wt,wt_test,
                                   alpha = 0.1,
                                   wthigh = 20, wtlow = 0.05,type,side){

    avg_wt <- mean(c(wt, wt_test))
    wt <- censoring(wt / avg_wt, wthigh, wtlow)
    wt_test <- censoring(wt_test / avg_wt, wthigh, wtlow)

    totw <- sum(wt)
    wt <- wt / totw
    qt <- (1 + wt_test / totw) * (1 - alpha)
    qt <- pmin(qt, 1)
    Yslack <- weightedConformalCutoff(Yscore, wt, qt)

    if (type == "CQR" && side == "two"){
        Ylo <- Yhat_test[, 1] - Yslack
        Yup <- Yhat_test[, 2] + Yslack
    } else if (type == "mean" && side == "two"){
        Ylo <- Yhat_test - Yslack
        Yup <- Yhat_test + Yslack
    } else if (side == "above"){
        Ylo <- -Inf
        Yup <- Yhat_test + Yslack
    } else if (side == "below"){
        Ylo <- Yhat_test - Yslack
        Yup <- Inf
    }

    res <- data.frame(lower = Ylo, upper = Yup)
    return(res)
}

guessClass <- function(x){
    if (length(unique(x)) == 2 || is.logical(x) || (is.factor(x) && nlevels(x) == 2)){
        dist <- "bernoulli"
    } else if (is.factor(x) && nlevels(x) > 2){
        dist <- "multinomial"
    } else if (is.numeric(x)){
        dist <- "gaussian"
    }
    return(dist)
}

quantRF <- function(Y, X, Xtest, quantiles, ...){
    fit <- grf::quantile_forest(X, Y, quantiles = quantiles, ...)
    res <- predict(fit, Xtest, quantiles = quantiles)
    res <- res[[1]]
    if (length(quantiles) == 1){
        res <- as.numeric(res)
    } else {
        res <- as.matrix(res)
    }
    return(res)
}

RF <- function(Y, X, Xtest, ...){
    dist <- guessClass(Y)
    if (dist == "gaussian"){
        fit <- randomForest::randomForest(x = X, y = Y, ...)
        res <- predict(fit, newdata = Xtest)
        res <- as.numeric(res)
    } else if (dist == "bernoulli"){
        if (!is.factor(Y)){
            Y <- as.factor(Y)
        }
        fit <- randomForest::randomForest(x = X, y = Y, ...)
        res <- predict(fit, newdata = Xtest, type = "prob")
        res <- as.numeric(res[, 2])
    } else if (dist == "multinomial"){
        if (!is.factor(Y)){
            Y <- as.factor(Y)
        }
        fit <- randomForest::randomForest(x = X, y = Y, ...)
        res <- predict(fit, newdata = Xtest, type = "prob")
        res <- as.matrix(res)
    }
    return(res)
}

quantBoosting <- function(Y, X, Xtest, quantiles, n.trees = 100, ...){
    # if (class(X) != "data.frame"){
    #     X <- as.data.frame(X)
    #     Xtest <- as.data.frame(Xtest)
    #     names(Xtest) <- names(X)
    # }
        X <- as.data.frame(X)
        Xtest <- as.data.frame(Xtest)
        names(Xtest) <- names(X)
    data <- data.frame(Y = Y, X)
    fit <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantiles[1]), data = data, n.trees = n.trees, ...)
    res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
    if (length(quantiles) == 2){
        fit2 <- gbm::gbm(Y ~ ., distribution = list(name = "quantile", alpha = quantiles[2]), data = data, n.trees = n.trees, ...)
        res2 <- predict(fit2, Xtest, type = "response", n.trees = n.trees)
        res <- cbind(res, res2)
    }
    return(res)
}

Boosting <- function(Y, X, Xtest, n.trees = 100, ...){
    # if (class(X) != "data.frame"){
    #     X <- as.data.frame(X)
    #     Xtest <- as.data.frame(Xtest)
    #     names(Xtest) <- names(X)
    # }
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
    distribution <- guessClass(Y)
    if (distribution == "bernoulli" && is.factor(Y)){
        Y <- as.numeric(Y) - 1
    }
    data <- data.frame(Y = Y, X)
    fit <- gbm::gbm(Y ~ ., distribution = distribution, data = data, n.trees = n.trees, ...)
    res <- predict(fit, Xtest, type = "response", n.trees = n.trees)
    if (distribution == "multinomial"){
        res <- matrix(res, nrow = nrow(Xtest))
    }
    return(res)
}

Logistic <- function(Y, X, Xtest, ...){
    # if (class(X) != "data.frame"){
    #     X <- as.data.frame(X)
    #     Xtest <- as.data.frame(Xtest)
    #     names(Xtest) <- names(X)
    # }
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
    distribution <- guessClass(Y)
    if (distribution == "bernoulli" && is.factor(Y)){
        Y <- as.numeric(Y) - 1
    }
    data <- data.frame(Y = Y, X)
    fit <- stats::glm(Y ~ ., family=binomial(link='logit'),data=data)
    res <- predict(fit, Xtest, type = "response")
    if (distribution == "multinomial"){
        res <- matrix(res, nrow = nrow(Xtest))
    }
    return(res)
}

PenLogistic <- function(Y, X, Xtest, ...){
    # if (class(X) != "data.frame"){
    #     X <- as.data.frame(X)
    #     Xtest <- as.data.frame(Xtest)
    #     names(Xtest) <- names(X)
    # }
    X <- as.matrix(X)
    Xtest <- as.matrix(Xtest)
    names(Xtest) <- names(X)
    distribution <- guessClass(Y)
    if (distribution == "bernoulli" && is.factor(Y)){
        Y <- as.numeric(Y) - 1
    }
    data <- data.frame(Y = Y, X)
    fit <- glmnet::cv.glmnet(X, Y, family="binomial",alpha = 1)
    res <- predict(fit, newx=Xtest, type = "response",s="lambda.min")
    if (distribution == "multinomial"){
        res <- matrix(res, nrow = nrow(Xtest))
    }
    return(res)
}

quantBART <- function(Y, X, Xtest, quantiles,
                      ndpost = 100, ...){
    # if (class(X) != "data.frame"){
    #     X <- as.data.frame(X)
    #     Xtest <- as.data.frame(Xtest)
    #     names(Xtest) <- names(X)
    # }
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
    fit <- bartMachine::bartMachine(X, Y, verbose = FALSE)
    if (length(quantiles) == 2){
        if (sum(quantiles) != 1){
            warning("Two quantiles should sum up to 1.")
        }
        ci_conf <- quantiles[2] - quantiles[1]
        res <- bartMachine::calc_prediction_intervals(
                                fit, new_data = Xtest,
                                pi_conf = 0.95)$interval
        res <- as.matrix(res)
    } else if (length(quantiles) == 1){
        if (quantiles[1] > 0.5){
            ci_conf <- 2 * quantiles[1]
            res <- bartMachine::calc_prediction_intervals(
                                    fit, new_data = Xtest,
                                    pi_conf = 0.95)$interval[, 2]
            res <- as.numeric(res)
        }  else{
            ci_conf <- 2 * (1 - quantiles[1])
            res <- bartMachine::calc_prediction_intervals(
                                    fit, new_data = Xtest,
                                    pi_conf = 0.95)$interval[, 1]
            res <- as.numeric(res)
        }
    }
    return(res)
}

BART <- function(Y, X, Xtest,
                 ndpost = 100, ...){
    # if (class(X) != "data.frame"){
    #     X <- as.data.frame(X)
    #     Xtest <- as.data.frame(Xtest)
    #     names(Xtest) <- names(X)
    # }
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
    fit <- bartMachine::bartMachine(X, Y, verbose = FALSE)
    res <- predict(fit, Xtest)
    return(res)
}

str_psfun <- function(method){
  if (method == "RF"){
    if (!requireNamespace("randomForest")){
      stop("randomForest package needs to be installed")
    }
    return(RF)
  } else if (method == "Boosting"){
    if (!requireNamespace("gbm")){
      stop("gbm package needs to be installed")
    }
    return(Boosting)
  } else {
    stop(paste0(method, " is not supported. Please input a valid string or a function that meets the minimal requirements described in the man page"))
  }
}

PropScore <- function(Ttrain,Xtrain,psfun, psparams){

  if (is.null(psfun)){
    psfun <- Boosting
  } else if (is.character(psfun)){
    psfun <- str_psfun(psfun[1])
  }

  psparams0 <- psparams
  psparams <- c(list(Y = Ttrain, X = Xtrain), psparams0)
  PropScorefun <- function(X){
    ps <- do.call(psfun, c(list(Xtest = X), psparams))
    ps
  }
  PropScorefun
}

str_outfun <- function(method){
  if (method == "RF"){
    if (!requireNamespace("randomForest")){
      stop("randomForest package needs to be installed")
    }
    return(RF)
  } else if (method == "quantRF"){
    if (!requireNamespace("grf")){
      stop("grf package needs to be installed")
    }
    return(quantRF)
  } else if (method == "Boosting"){
    if (!requireNamespace("gbm")){
      stop("gbm package needs to be installed")
    }
    return(Boosting)
  } else if (method == "quantBoosting"){
    if (!requireNamespace("gbm")){
      stop("gbm package needs to be installed")
    }
    return(quantBoosting)
  } else if (method == "BART"){
    if (!requireNamespace("bartMachine")){
      stop("bartMachine package needs to be installed")
    }
    return(BART)
  } else if (method == "quantBART"){
    if (!requireNamespace("bartMachine")){
      stop("bartMachine package needs to be installed")
    }
    return(quantBART)
  } else {
    stop(paste0(method, " is not supported. Please input a valid string or a function that meets the minimal requirements described in the man page"))
  }
}

QuantileReg <- function(Ytrain,Xtrain,outfun,outparams=NULL,quantiles){
  outparams <- c(outparams, list(quantiles = quantiles))
  outparams <- c(list(Y = Ytrain, X = Xtrain), outparams)

  Ymodel <- function(X){
    fit <- do.call(outfun, c(outparams, list(Xtest = X)))
    fit
  }
  Ymodel
}

SampleSplit <- function(Y,trainprop){
  T <- as.numeric(!is.na(Y))
  inds1 <- which(T == 1)
  inds0 <- which(T == 0)
  n1 <- length(inds1)
  n0 <- length(inds0)
  trainid1 <- sample(n1, floor(n1 * trainprop))
  trainid0 <- sample(n0, floor(n0 * trainprop))
  trainid <- c(inds1[trainid1], inds0[trainid0])
  return(list(trainid=trainid, trainid1=trainid1,trainid0=trainid0))
}

Naive <- function(X,Y,T,Deta,Xtest,trainprop,psfun,psparams,outfun,outparams,quantiles,type,side,alpha,wthigh,wtlow){

  n <- length(Y)
  #split samples for training and calibration
  YY <- Y
  YY[!T] <- NA
  T <- as.numeric(!is.na(YY))
  res <- SampleSplit(YY,trainprop)
  #res$trainid: index for all training
  #res$trainid1: index for all training with T=1
  #!!!need to select T=1 first and then select trainid1
  #res$trainid0: index for all training with T=0

  #estimate the propensity score
  Xtrain <- X[res$trainid, ,drop=FALSE]
  Ttrain <- T[res$trainid]
  PSFit_T <- PropScore(Ttrain,Xtrain,psfun, psparams)

  #estimate the propensity score for binary A
  tp_T <- PSFit_T(X)
  tp_D <- exp(cbind(1,X)%*%Deta)/(1+exp(cbind(1,X)%*%Deta))
  ps_T_test <- PSFit_T(Xtest)

  #random sample from binary varibale(D)
  D <- rbinom(n = n, size = 1, prob=tp_D)
  O <- I(T==D)
  Nts <- sum(O)
  #quantile regression fitting after selecting A=T
  X1 <- X[res$trainid, ,drop=FALSE]
  Y1 <- Y[res$trainid]
  #select O=1 from training data
  ind <- which(O[res$trainid] == 1)
  Xtrain <- X1[ind, ,drop=FALSE]
  Ytrain <- Y1[ind]
  OutputFit <- QuantileReg(Ytrain,Xtrain,outfun,outparams,quantiles)
  Yhat_test <- OutputFit(Xtest)

  #data for validation
  X0 <- X[-res$trainid, ,drop=FALSE]
  Y0 <- Y[-res$trainid]
  ind <- which(O[-res$trainid] == 1)
  Xval <- X0[ind, ,drop=FALSE]
  Yval <- Y0[ind]
  ps_T_val <- PSFit_T(as.data.frame(Xval))
  Yhat_val <- OutputFit(Xval)

  Yscore <- conformalScore(Yval, Yhat_val, type, side)
  #construct the confidence interval by weighted CP
  wt <- rep(1,length(ps_T_val))
  wt_test <- rep(1,length(ps_T_test))
  PredInt<- predict.conformalSplit(Yhat_test, Yscore,wt,wt_test,alpha,wthigh,wtlow,type,side)

  lefts <- PredInt[,1]
  rights <- PredInt[,2]
  Lens <- rights-lefts

  return(list(lefts=lefts, rights=rights, Lens=Lens, Nts=Nts))
}

