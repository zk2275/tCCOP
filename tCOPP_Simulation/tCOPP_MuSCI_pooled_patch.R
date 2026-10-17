## =====================================================================
## tCOPP_MuSCI_pooled_patch.R  (AUTO-GENERATED -- do not hand-edit)
## ---------------------------------------------------------------------
## Overrides MuSCI_ConformalCI_MuS() and MuSCI_ConformalCI() so that the
## POOLED estimator [tCOPP-EIF(P)] is returned in addition to the federated
## / target estimator. Source AFTER tCOPP_functions.R:
##     source("tCOPP_functions.R")
##     source("tCOPP_MuSCI_pooled_patch.R")
##
## Backward compatible: out$lefts/out$rights/out$Lens are unchanged
##   (federated for _MuS, target-site for single-source). The pooled
##   interval is exposed as out$lefts_pool/out$rights_pool/out$Lens_pool.
## =====================================================================

## Majority-vote (B>1) / single-split (B==1) aggregation, byte-identical to the
## loop in the original MuSCI functions, but operating on a generic pair of
## per-split endpoint matrices (rows = splits, cols = test points).
.tCOPP_aggregate <- function(YBleft, YBright, B) {
  lefts <- NULL; rights <- NULL; Lens <- NULL
  if (B > 1) {
    ntest <- ncol(YBleft)
    for (j in 1:ntest) {
      Ys <- c(YBleft[, j], YBright[, j])
      Ws <- c(rep(1, nrow(YBleft)), rep(1, nrow(YBright)))
      Ss <- c(rep(1, nrow(YBleft)), rep(0, nrow(YBright)))
      tol <- B / 2
      ind <- order(Ys); Ys <- Ys[ind]; Ws <- Ws[ind]; Ss <- Ss[ind]
      count <- 0; leftpoints <- NULL; rightpoints <- NULL
      for (i in 1:length(Ys)) {
        if (Ss[i] == 1) {
          count <- count + Ws[i]
          if ((count > tol) & (count - Ws[i] <= tol)) leftpoints <- c(leftpoints, Ys[i])
        } else if (Ss[i] == 0) {
          if ((count > tol) & (count - Ws[i] <= tol)) rightpoints <- c(rightpoints, Ys[i])
          count <- count - Ws[i]
        }
      }
      lefts[[j]] <- leftpoints
      rights[[j]] <- rightpoints
      Lens[j] <- sum(rightpoints - leftpoints)
    }
  } else if (B == 1) {
    ntest <- ncol(YBleft)
    for (j in 1:ntest) {
      lefts[[j]]  <- YBleft[1, j]
      rights[[j]] <- YBright[1, j]
      Lens[j]     <- YBright[1, j] - YBleft[1, j]
    }
  }
  list(lefts = lefts, rights = rights, Lens = Lens)
}

## ---- MuSCI_ConformalCI_MuS  (federated default; pooled exposed) ----
MuSCI_ConformalCI_MuS <- function(X,Y,X_ex,Y_ex,T,T_ex,Deta,Xtest,method,B,trainprop,psfun,psparams,outfun,outparams,quantiles,type,side,weighted,
####                  conf.score="ASR",
                  conf.score,
                  ps.Library="SL.glm",
                  m.Library="SL.glm",
                  p1,#0.5
                  seed=10012) {
  YBleft <- NULL
  YBright <- NULL
  YBleft_fed  <- NULL
  YBright_fed <- NULL
  YBleft_pool <- NULL
  YBright_pool <- NULL
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
    phi[ind1] <- (m_s1_val[D2_s1$R==0]-(1-alpha))/mean(D2$Site==1 & D2$R==0)
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
    YBleft_fed   <- rbind(YBleft_fed,   PI.fed2$lower)
    YBright_fed  <- rbind(YBright_fed,  PI.fed2$upper)
    YBleft_pool  <- rbind(YBleft_pool,  PI.pool$lower)
    YBright_pool <- rbind(YBright_pool, PI.pool$upper)
  }

  # [pooled-patch] aggregate federated and pooled intervals separately
  agg_fed  <- .tCOPP_aggregate(YBleft_fed,  YBright_fed,  B)
  agg_pool <- .tCOPP_aggregate(YBleft_pool, YBright_pool, B)
  return(list(lefts = agg_fed$lefts, rights = agg_fed$rights, Lens = agg_fed$Lens,
              lefts_pool = agg_pool$lefts, rights_pool = agg_pool$rights, Lens_pool = agg_pool$Lens,
              Nts = Nts, r0 = r0))

}


## ---- MuSCI_ConformalCI  (target-site default; pooled exposed) ----
MuSCI_ConformalCI <- function(X,Y,T,Deta,Xtest,method,B,trainprop,psfun,psparams,outfun,outparams,quantiles,type,side,weighted,
####                  conf.score="ASR",
                  conf.score,
                  ps.Library="SL.glm",
                  m.Library="SL.glm",
                  p1,#0.5
                  seed=10012) {
  YBleft <- NULL
  YBright <- NULL
  YBleft_tgt  <- NULL
  YBright_tgt <- NULL
  YBleft_pool <- NULL
  YBright_pool <- NULL
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
    PI.pool <- data.frame(lwr= -rhat_pool+pred_y0, upr=rhat_pool+pred_y0)
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
    PI.pool <- data.frame(lwr= -rhat_pool*rho0+pred_y0, upr=rhat_pool*rho0+pred_y0)

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
    PI.pool <- data.frame(lwr= pred_y0.lw-rhat_pool, upr= pred_y0.up+rhat_pool)

  }
  colnames(PI.tgSt) <- c("lower", "upper")
  colnames(PI.pool) <- c("lower", "upper")
    #PredInt[,1]=PI.tgSt$lower
    #PredInt[,2]=PI.tgSt$upper
    #YBleft <- rbind(YBleft,PredInt[,1])
    #YBright <- rbind(YBright,PredInt[,2])
    YBleft_tgt   <- rbind(YBleft_tgt,   PI.tgSt$lower)
    YBright_tgt  <- rbind(YBright_tgt,  PI.tgSt$upper)
    YBleft_pool  <- rbind(YBleft_pool,  PI.pool$lower)
    YBright_pool <- rbind(YBright_pool, PI.pool$upper)
  }

  # [pooled-patch] aggregate target-site and pooled intervals separately
  agg_tgt  <- .tCOPP_aggregate(YBleft_tgt,  YBright_tgt,  B)
  agg_pool <- .tCOPP_aggregate(YBleft_pool, YBright_pool, B)
  return(list(lefts = agg_tgt$lefts, rights = agg_tgt$rights, Lens = agg_tgt$Lens,
              lefts_pool = agg_pool$lefts, rights_pool = agg_pool$rights, Lens_pool = agg_pool$Lens,
              Nts = Nts, r0 = r0))

}

