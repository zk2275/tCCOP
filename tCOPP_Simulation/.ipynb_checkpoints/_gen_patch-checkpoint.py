#!/usr/bin/env python3
"""
Generate tCOPP_MuSCI_pooled_patch.R.

Strategy: copy the two MuSCI function bodies VERBATIM from tCOPP_functions.R
(so 99% of the code is byte-identical to the tested original) and apply only
a handful of surgical, slice-scoped edits so that:

  * MuSCI_ConformalCI_MuS  now returns the POOLED interval (PI.pool) alongside
    the federated interval (PI.fed2). Federated stays the default (backward
    compatible), pooled is exposed as *_pool.
  * MuSCI_ConformalCI (single source) now returns the POOLED interval
    (built from rhat_pool) alongside the target-site interval (PI.tgSt).
    Target stays the default (backward compatible), pooled exposed as *_pool.

The interval-aggregation logic (majority vote for B>1, single split for B==1)
is factored into one helper, .tCOPP_aggregate(), that reproduces the original
loop exactly and is applied once per interval.

Sourcing this file AFTER tCOPP_functions.R overrides the two functions
(R last-definition-wins).
"""

SRC = "tCOPP_functions.R"
OUT = "tCOPP_MuSCI_pooled_patch.R"

with open(SRC) as f:
    lines = f.readlines()

def slice_lines(a, b):
    # 1-indexed inclusive
    return "".join(lines[a-1:b])

mus  = slice_lines(36, 687)    # MuSCI_ConformalCI_MuS  <- function(...) { ... }
sing = slice_lines(737, 1112)  # MuSCI_ConformalCI      <- function(...) { ... }

# ---------------------------------------------------------------------------
# helper reproducing the ORIGINAL aggregation, operating on generic YB matrices
# ---------------------------------------------------------------------------
helper = r'''## =====================================================================
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

'''

# ---------------------------------------------------------------------------
# EDIT 1: MuSCI_ConformalCI_MuS  (federated + pooled)
# ---------------------------------------------------------------------------
# 1a. add pooled accumulators next to the fed accumulators (init block)
mus = mus.replace(
    "  YBleft <- NULL\n  YBright <- NULL\n  Nts <- NULL\n  r0<-NULL\n  r1<-NULL\n",
    "  YBleft <- NULL\n  YBright <- NULL\n"
    "  YBleft_fed  <- NULL\n  YBright_fed <- NULL\n"
    "  YBleft_pool <- NULL\n  YBright_pool <- NULL\n"
    "  Nts <- NULL\n  r0<-NULL\n  r1<-NULL\n",
    1)

# 1b. accumulate BOTH fed2 and pool inside the b-loop
mus = mus.replace(
    "    YBleft <- rbind(YBleft,PI.fed2$lower)\n"
    "    YBright <- rbind(YBright,PI.fed2$upper)\n",
    "    YBleft_fed   <- rbind(YBleft_fed,   PI.fed2$lower)\n"
    "    YBright_fed  <- rbind(YBright_fed,  PI.fed2$upper)\n"
    "    YBleft_pool  <- rbind(YBleft_pool,  PI.pool$lower)\n"
    "    YBright_pool <- rbind(YBright_pool, PI.pool$upper)\n",
    1)

# 1c. replace the whole post-loop tail (aggregation + return) with helper calls
mus_tail_old = slice_lines(633, 685)   # exact original tail text (within _MuS)
mus_tail_new = (
    "  # [pooled-patch] aggregate federated and pooled intervals separately\n"
    "  agg_fed  <- .tCOPP_aggregate(YBleft_fed,  YBright_fed,  B)\n"
    "  agg_pool <- .tCOPP_aggregate(YBleft_pool, YBright_pool, B)\n"
    "  return(list(lefts = agg_fed$lefts, rights = agg_fed$rights, Lens = agg_fed$Lens,\n"
    "              lefts_pool = agg_pool$lefts, rights_pool = agg_pool$rights, Lens_pool = agg_pool$Lens,\n"
    "              Nts = Nts, r0 = r0))\n"
)
assert mus_tail_old in mus, "MuS tail anchor not found"
mus = mus.replace(mus_tail_old, mus_tail_new, 1)

# ---------------------------------------------------------------------------
# EDIT 2: MuSCI_ConformalCI  (target-site + pooled)
# ---------------------------------------------------------------------------
# 2a. add pooled accumulators (init block)
sing = sing.replace(
    "  YBleft <- NULL\n  YBright <- NULL\n  Nts <- NULL\n  r0<-NULL\n  r1<-NULL\n",
    "  YBleft <- NULL\n  YBright <- NULL\n"
    "  YBleft_tgt  <- NULL\n  YBright_tgt <- NULL\n"
    "  YBleft_pool <- NULL\n  YBright_pool <- NULL\n"
    "  Nts <- NULL\n  r0<-NULL\n  r1<-NULL\n",
    1)

# 2b. build PI.pool in each of the three conformal-score branches (from rhat_pool)
sing = sing.replace(
    "    PI.tgSt <- data.frame(lwr= -rhat_tgSt+pred_y0, upr=rhat_tgSt+pred_y0)\n",
    "    PI.tgSt <- data.frame(lwr= -rhat_tgSt+pred_y0, upr=rhat_tgSt+pred_y0)\n"
    "    PI.pool <- data.frame(lwr= -rhat_pool+pred_y0, upr=rhat_pool+pred_y0)\n",
    1)
sing = sing.replace(
    "    PI.tgSt <- data.frame(lwr= -rhat_tgSt*rho0+pred_y0, upr=rhat_tgSt*rho0+pred_y0)\n",
    "    PI.tgSt <- data.frame(lwr= -rhat_tgSt*rho0+pred_y0, upr=rhat_tgSt*rho0+pred_y0)\n"
    "    PI.pool <- data.frame(lwr= -rhat_pool*rho0+pred_y0, upr=rhat_pool*rho0+pred_y0)\n",
    1)
sing = sing.replace(
    "    PI.tgSt <- data.frame(lwr= pred_y0.lw-rhat_tgSt, upr= pred_y0.up+rhat_tgSt)\n",
    "    PI.tgSt <- data.frame(lwr= pred_y0.lw-rhat_tgSt, upr= pred_y0.up+rhat_tgSt)\n"
    "    PI.pool <- data.frame(lwr= pred_y0.lw-rhat_pool, upr= pred_y0.up+rhat_pool)\n",
    1)

# 2c. rename columns of BOTH PI.tgSt and PI.pool
sing = sing.replace(
    '  colnames(PI.tgSt) <- c("lower", "upper")\n',
    '  colnames(PI.tgSt) <- c("lower", "upper")\n'
    '  colnames(PI.pool) <- c("lower", "upper")\n',
    1)

# 2d. accumulate BOTH target and pool inside the b-loop
sing = sing.replace(
    "    YBleft <- rbind(YBleft,PI.tgSt$lower)\n"
    "    YBright <- rbind(YBright,PI.tgSt$upper)\n",
    "    YBleft_tgt   <- rbind(YBleft_tgt,   PI.tgSt$lower)\n"
    "    YBright_tgt  <- rbind(YBright_tgt,  PI.tgSt$upper)\n"
    "    YBleft_pool  <- rbind(YBleft_pool,  PI.pool$lower)\n"
    "    YBright_pool <- rbind(YBright_pool, PI.pool$upper)\n",
    1)

# 2e. replace the post-loop tail (aggregation + return) with helper calls
sing_tail_old = slice_lines(1060, 1110)  # exact original tail text (within single)
sing_tail_new = (
    "  # [pooled-patch] aggregate target-site and pooled intervals separately\n"
    "  agg_tgt  <- .tCOPP_aggregate(YBleft_tgt,  YBright_tgt,  B)\n"
    "  agg_pool <- .tCOPP_aggregate(YBleft_pool, YBright_pool, B)\n"
    "  return(list(lefts = agg_tgt$lefts, rights = agg_tgt$rights, Lens = agg_tgt$Lens,\n"
    "              lefts_pool = agg_pool$lefts, rights_pool = agg_pool$rights, Lens_pool = agg_pool$Lens,\n"
    "              Nts = Nts, r0 = r0))\n"
)
assert sing_tail_old in sing, "single-source tail anchor not found"
sing = sing.replace(sing_tail_old, sing_tail_new, 1)

with open(OUT, "w") as f:
    f.write(helper)
    f.write("## ---- MuSCI_ConformalCI_MuS  (federated default; pooled exposed) ----\n")
    f.write(mus)
    f.write("\n\n")
    f.write("## ---- MuSCI_ConformalCI  (target-site default; pooled exposed) ----\n")
    f.write(sing)
    f.write("\n")

print("Wrote", OUT)
