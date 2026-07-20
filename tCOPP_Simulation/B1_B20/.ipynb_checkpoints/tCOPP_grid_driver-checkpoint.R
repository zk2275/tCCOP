## =====================================================================
## tCOPP_grid_driver.R
## ---------------------------------------------------------------------
## Defines the full simulation grid for the 60-boxplot deliverable:
##
##   3 scenarios  x  5 methods  x  2 source-modes  x  2 values of B  = 60 cells
##
##   scenarios : (a) CCOD holds for ALL sources
##               (b) CCOD violated for external sources (simple)
##               (c) CCOD violated for external sources (complex, paper Table)
##   methods   : COPP, tCOPP-W, tCOPP-EIF(I), tCOPP-EIF(M), tCOPP-EIF(P)
##   modes     : "internal"  (internal source D1 only)
##               "multi"     (D1 + external D2..DK)
##   B         : 1  and  20   (multi-sampling / majority-vote splits, B_MS)
##
## Source order matters (last definition wins):
##     source("tCOPP_functions.R")            # shared library
##     source("tCOPP_MuSCI_pooled_patch.R")   # exposes pooled EIF(P) intervals
##     source("tCOPP_grid_driver.R")          # this file
##
## The data-generating code and the CP/AL computation are copied from the
## existing figure files (Fig_2a/2b Multi-Source) so numbers stay comparable
## to the current manuscript figures. Only the external-site outcome model
## (genY_singlestage_ex_*) changes across scenarios (a)/(b)/(c).
##
## Interpretation of the two source-modes for each method
## ------------------------------------------------------
##   COPP / tCOPP-W / tCOPP-EIF(I):
##       internal -> fit on D1 only
##       multi    -> fit on the NAIVE POOL of D1..DK stacked as one source
##                   (this is the comparison that is expected to break under
##                    CCOD violation -- exactly the point of scenarios b,c)
##   tCOPP-EIF(M) (federated):
##       internal -> reduces mathematically to tCOPP-EIF(I) (single site)
##       multi    -> federated aggregation over D1 (internal) + D2..DK
##   tCOPP-EIF(P) (pooled oracle, site indicator as covariate):
##       internal -> pooled estimator on {D1, D0}          (single-source rhat_pool)
##       multi    -> pooled estimator on {D1..DK, D0}       (MuS rhat_pool)
##
## Efficiency: within one replicate we fit each model ONCE and read off every
## interval it produces (e.g. the multi-source EIF fit yields BOTH the
## federated and the pooled interval), so paired cells share the same data and
## the same fit -- consistent with the project's scientific-integrity rule.
## =====================================================================

## ---------------------------------------------------------------------
## Global simulation parameters (match the existing figure files)
## ---------------------------------------------------------------------
GRID_PARAMS <- list(
  d          = 4,
  trainprop  = 0.75,
  quantiles  = c(0.05, 0.95),
  type       = "CQR",     # conformal-score type for ConformalCI (COPP / tCOPP-W)
  side       = "two",
  alpha      = 0.1,       # nominal 1-alpha = 0.9
  wthigh     = 20,
  wtlow      = 0.05,
  conf.score = "localASR",# conformal score for the EIF methods
  Ntl        = 15000,     # super-population size before RCT selection
  ntest      = 10000,     # target-sample size
  ex_site_num  = 4,       # number of EXTERNAL sources (K-1); total sources K = 5
  ex_site_size = 2000,    # per external source (paper uses ~2000; Fig2b-pooled used 500)
  w.method_tw  = 2,       # covariate-shift weighting for tCOPP-W (>=2 = real weights)
  w.method_copp= 1,       # COPP uses unit weights (no covariate-shift correction)
  Teta   = c(-0.5, -0.5, -0.5, -0.5, -0.5),          # behavior-policy logits
  Deta   = c(-0.5,  1,    1,   -1,   -1),            # policy-shift logits (pi_e vs pi_b)
  alphas = as.matrix(c(1.35, -4, -4, 0, 0), 5, 1),   # RCT selection: sigmoid(1.35-4X1-4X2)
  outfun = NULL,   # set below to quantRF (needs functions sourced first)
  psfun  = NULL,   # set below to Logistic
  seed_internal = 10012    # passed straight through to the MuSCI functions
)

## psfun / outfun require the shared library to be sourced already.
if (exists("Logistic")) GRID_PARAMS$psfun  <- Logistic
if (exists("quantRF"))  GRID_PARAMS$outfun <- quantRF

## IMPORTANT: MuSCI_ConformalCI() and MuSCI_ConformalCI_MuS() read `alpha` as a
## FREE variable via lexical scope -- it is NOT one of their formal arguments.
## The original figure chunks satisfy this with a top-level `alpha <- 0.1`.
## We do the same here so the estimating-equation root-finder inside those
## functions can resolve `alpha` in every (forked) worker. Without this the
## EIF cells silently return NA ("object 'alpha' not found").
alpha <- GRID_PARAMS$alpha

## The five methods, in left-to-right plotting order, each shown in its NATURAL
## source-mode (mixed-axis layout matching the manuscript Figure 2 caption
## "internal source only and multiple sources"):
##   COPP / tCOPP-W / tCOPP-EIF(I)  -> internal source D1 only
##   tCOPP-EIF(M) / tCOPP-EIF(P)    -> multiple sources (D1 + external)
GRID_METHODS <- c("COPP", "tCOPP-W", "tCOPP-EIF(I)", "tCOPP-EIF(M)", "tCOPP-EIF(P)")
GRID_CELLS   <- data.frame(
  method = c("COPP", "tCOPP-W", "tCOPP-EIF(I)", "tCOPP-EIF(M)", "tCOPP-EIF(P)"),
  mode   = c("internal", "internal", "internal", "multi", "multi"),
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------------
## Scenario-specific EXTERNAL-site outcome models
##   internal source D1 and target D0 always use genY_singlestage() (CCOD holds
##   internal<->target by construction). Only the external sites differ.
## ---------------------------------------------------------------------

# (a) CCOD holds for ALL sources: external sites share the internal model
#     (coefficient 1 on X3^3, identical to genY_singlestage) -> no shift.
genY_ex_a <- function(X, T, site_num) {
  n <- length(T)
  1 + X[,1] - X[,2] + X[,3]^{3} + exp(X[,4]) +
    T*(3 - 5*X[,1] + 2*X[,2] - 3*X[,3] + X[,4]) +
    (1 + T)*(1 + X[,1] + X[,2] + X[,3] + X[,4])*rnorm(n)
}

# (b) CCOD violated (simple): external site k scales the X3^3 term by site_num,
#     so its conditional outcome law differs from the target (coeff k != 1).
genY_ex_b <- function(X, T, site_num) {
  n <- length(T)
  1 + X[,1] - X[,2] + site_num*X[,3]^{3} + exp(X[,4]) +
    T*(3 - 5*X[,1] + 2*X[,2] - 3*X[,3] + X[,4]) +
    (1 + T)*(1 + X[,1] + X[,2] + X[,3] + X[,4])*rnorm(n)
}

# (c) CCOD violated (complex): site-specific f_k, g_k, h_k per manuscript Table 4.
#     site_num is the site index k in {2,3,4,5}.
genY_ex_c <- function(X, T, site_num) {
  n  <- length(T)
  X1 <- X[,1]; X2 <- X[,2]; X3 <- X[,3]; X4 <- X[,4]
  k  <- site_num
  if (k == 2) {
    f <-  2*X2 + X3^2 + 4*exp(X4)
    g <-  3 - 5*X1 - 3*X3 + 4*X4
    h <-  1 - 2*X1 + X2 + X3 + 4*X4
  } else if (k == 3) {
    f <-  2*X2 + X3^2 - 4*exp(X4)
    g <-  3 + 5*X1 + 2*X2 + 3*X3 - 4*X4
    h <-  1 + X1 + X2 - 3*X3 + X4
  } else if (k == 4) {
    f <- -X2 - 3*X3^2 + 4*exp(X4)
    g <-  3 + 5*X1 + 2*X2 + 3*X3 + 4*X4
    h <-  1 + X1 + X2 - 3*X3 + 4*X4
  } else {  # k == 5 (and any k>5 fallback)
    f <-  2*X2 - 3*X3^2 - 4*exp(X4)
    g <-  3 - 5*X1 - 3*X3 - 4*X4
    h <-  1 - 2*X1 + X2 + X3 + X4
  }
  1 + X1 + f + T*g + (1 + T)*h*rnorm(n)
}

scenario_genY_ex <- function(scenario) {
  switch(scenario,
         a = genY_ex_a,
         b = genY_ex_b,
         c = genY_ex_c,
         stop("scenario must be 'a', 'b', or 'c'"))
}

## ---------------------------------------------------------------------
## CP / AL for one fitted 'out' object, faithful to the figure-file logic.
##   which = "default" -> out$lefts / out$rights / out$Lens
##   which = "pool"    -> out$lefts_pool / out$rights_pool / out$Lens_pool
## The B>1 branch reproduces the exact convention used in the existing files.
## ---------------------------------------------------------------------
cp_al <- function(out, which, Ytest, ntest, B) {
  if (which == "pool") {
    lefts <- out$lefts_pool; rights <- out$rights_pool; Lens <- out$Lens_pool
  } else {
    lefts <- out$lefts;      rights <- out$rights;      Lens <- out$Lens
  }
  if (B == 1) {
    L <- unlist(lefts); R <- unlist(rights)
    ind <- (Ytest < R) & (Ytest > L)
    CP <- sum(ind) / ntest
    AL <- mean(unlist(Lens))
  } else {
    ind <- logical(ntest)
    for (kk in 1:ntest) {
      NInt <- length(lefts[kk])   # matches existing files' convention (single bracket)
      flag <- unlist(lapply(seq_len(NInt), function(x)
        (lefts[[kk]][x] < Ytest[kk]) & (Ytest[kk] < rights[[kk]][x])))
      ind[kk] <- (sum(flag) > 0)
    }
    CP <- sum(ind) / ntest
    AL <- mean(unlist(Lens))
  }
  c(CP = CP, AL = AL)
}

## ---------------------------------------------------------------------
## One Monte-Carlo replicate for a given scenario & B.
## Returns a data.frame with 10 rows: one per (method, mode) cell.
## ---------------------------------------------------------------------
run_one_rep_grid <- function(rep_id, scenario, B, P = GRID_PARAMS,
                             base_seed = 20260701L) {

  ## reproducible, distinct RNG stream per (scenario, B, rep)
  set.seed(base_seed + 1000L*match(scenario, c("a","b","c")) + 7L*B + rep_id)

  d          <- P$d
  ntest      <- P$ntest
  Ntl        <- P$Ntl
  trainprop  <- P$trainprop
  quantiles  <- P$quantiles
  type       <- P$type
  side       <- P$side
  alpha      <- P$alpha
  wthigh     <- P$wthigh
  wtlow      <- P$wtlow
  conf.score <- P$conf.score
  psfun      <- P$psfun
  outfun     <- P$outfun
  Teta       <- P$Teta
  Deta       <- P$Deta
  alphas     <- P$alphas
  method     <- "BootSweight"           # subsampling estimator used throughout
  genY_ex    <- scenario_genY_ex(scenario)

  ## ---- super-population + internal RCT (verbatim structure) ----------
  X_N   <- matrix(runif(Ntl * d), nrow = Ntl, ncol = d)
  tp_T_N<- exp(cbind(1, X_N) %*% Teta) / (1 + exp(cbind(1, X_N) %*% Teta))
  T_N   <- rbinom(Ntl, 1, tp_T_N)
  Y_N   <- genY_singlestage(X_N, T_N)

  ## target sample D0 (constant selection; covariate law = U(0,1))
  Xtest    <- matrix(runif(ntest * d), nrow = ntest, ncol = d)
  tp_D_test<- exp(cbind(1, Xtest) %*% Deta) / (1 + exp(cbind(1, Xtest) %*% Deta))
  Dtest    <- rbinom(ntest, 1, tp_D_test)
  Ytest    <- genY_singlestage(Xtest, Dtest)

  ## external sources D2..DK (site index i+1); scenario-specific outcomes
  X_ex <- NULL; T_ex <- NULL; Y_ex <- NULL
  for (i in 1:P$ex_site_num) {
    Xs  <- matrix(runif(P$ex_site_size * d), nrow = P$ex_site_size, ncol = d)
    Xdf <- data.frame(Xs); Xdf$Site <- i + 1
    tpT <- exp(cbind(1, Xs) %*% Teta) / (1 + exp(cbind(1, Xs) %*% Teta))
    Ts  <- rbinom(P$ex_site_size, 1, tpT)
    Ys  <- genY_ex(Xs, Ts, i + 1)
    if (i == 1) { X_ex <- Xdf; T_ex <- Ts; Y_ex <- Ys }
    else        { X_ex <- rbind(X_ex, Xdf); T_ex <- c(T_ex, Ts); Y_ex <- c(Y_ex, Ys) }
  }

  ## internal source D1 via RCT selection sigmoid(1.35 - 4X1 - 4X2)
  delta      <- rbinom(Ntl, 1, delta.tilde.prob(cbind(1, X_N), alphas))
  sel        <- which(delta == 1)
  X          <- X_N[sel, ]
  T          <- T_N[sel]
  Y          <- Y_N[sel]

  Xtest_df <- as.data.frame(Xtest)

  ## helper: on a single cell error, return NA rather than aborting the whole rep
  safe <- function(expr) tryCatch(expr, error = function(e) {
    message("  cell error: ", conditionMessage(e)); NULL })

  ## named skeleton, one slot per (method.mode) cell
  rows <- setNames(vector("list", nrow(GRID_CELLS)),
                   paste(GRID_CELLS$method, GRID_CELLS$mode, sep = "."))

  ## ---- COPP (internal): weighted CP, unit covariate weights ----------
  rows[["COPP.internal"]] <- safe({
    wttl  <- learn_weights(X, Xtest, P$w.method_copp, Ntl)$w.tilde.rct
    wttlT <- learn_weights(X, Xtest, P$w.method_copp, Ntl)$w.tilde.rwe
    out <- ConformalCI(as.matrix(X), Y, T, Deta, Xtest_df, method, B, trainprop,
                       psfun, NULL, outfun, NULL, quantiles, type, side,
                       weighted = TRUE, alpha, wthigh, wtlow, wttl, wttlT)
    m <- cp_al(out, "default", Ytest, ntest, B)
    data.frame(method = "COPP", mode = "internal", CP = m["CP"], AL = m["AL"], row.names = NULL)
  })

  ## ---- tCOPP-W (internal): weighted CP, real covariate-shift weights -
  rows[["tCOPP-W.internal"]] <- safe({
    wttl  <- learn_weights(cbind(1, X), cbind(1, Xtest), P$w.method_tw, Ntl)$w.tilde.rct
    wttlT <- learn_weights(cbind(1, X), cbind(1, Xtest), P$w.method_tw, Ntl)$w.tilde.rwe
    out <- ConformalCI(as.matrix(X), Y, T, Deta, Xtest_df, method, B, trainprop,
                       psfun, NULL, outfun, NULL, quantiles, type, side,
                       weighted = TRUE, alpha, wthigh, wtlow, wttl, wttlT)
    m <- cp_al(out, "default", Ytest, ntest, B)
    data.frame(method = "tCOPP-W", mode = "internal", CP = m["CP"], AL = m["AL"], row.names = NULL)
  })

  ## ---- tCOPP-EIF(I) (internal): single-source EIF, target-site interval
  rows[["tCOPP-EIF(I).internal"]] <- safe({
    out <- MuSCI_ConformalCI(as.matrix(X), Y, T, Deta, Xtest_df, method, B, trainprop,
                             psfun, NULL, outfun, NULL, quantiles, type, side,
                             weighted = TRUE, conf.score,
                             ps.Library = "SL.glm", m.Library = "SL.glm",
                             p1 = trainprop, seed = P$seed_internal)
    m <- cp_al(out, "default", Ytest, ntest, B)
    data.frame(method = "tCOPP-EIF(I)", mode = "internal", CP = m["CP"], AL = m["AL"], row.names = NULL)
  })

  ## ---- tCOPP-EIF(M) & tCOPP-EIF(P) (multi): ONE multi-source EIF fit --
  ##      federated interval -> EIF(M);  pooled interval -> EIF(P)
  eif_mus <- safe(
    MuSCI_ConformalCI_MuS(as.matrix(X), Y, X_ex, Y_ex, T, T_ex, Deta, Xtest_df,
                          method, B, trainprop, psfun, NULL, outfun, NULL,
                          quantiles, type, side, weighted = TRUE, conf.score,
                          ps.Library = "SL.glm", m.Library = "SL.glm",
                          p1 = trainprop, seed = P$seed_internal)
  )
  if (!is.null(eif_mus)) {
    mM <- cp_al(eif_mus, "default", Ytest, ntest, B)   # federated = EIF(M)
    mP <- cp_al(eif_mus, "pool",    Ytest, ntest, B)   # pooled    = EIF(P)
    rows[["tCOPP-EIF(M).multi"]] <- data.frame(method = "tCOPP-EIF(M)", mode = "multi",
                                               CP = mM["CP"], AL = mM["AL"], row.names = NULL)
    rows[["tCOPP-EIF(P).multi"]] <- data.frame(method = "tCOPP-EIF(P)", mode = "multi",
                                               CP = mP["CP"], AL = mP["AL"], row.names = NULL)
  }

  ## assemble; fill any missing cell with NA so every rep returns 5 rows
  present <- rows[!vapply(rows, is.null, logical(1))]
  got <- if (length(present)) do.call(rbind, present) else
           data.frame(method = character(0), mode = character(0),
                      CP = numeric(0), AL = numeric(0))
  out_df <- merge(GRID_CELLS, got, by = c("method", "mode"), all.x = TRUE)
  out_df$scenario <- scenario
  out_df$B        <- B
  out_df$rep      <- rep_id
  out_df$method   <- factor(out_df$method, levels = GRID_METHODS)   # plotting order
  out_df <- out_df[order(out_df$method), ]
  out_df$method   <- as.character(out_df$method)
  out_df[, c("scenario", "B", "rep", "method", "mode", "CP", "AL")]
}
