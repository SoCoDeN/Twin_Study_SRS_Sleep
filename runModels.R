# ===========================
# Bivariate sex-limitation ACE (Quant then Full)
# Traits: rsrs_imp, rsds_imp
# Groups: MZf, DZf, MZm, DZm, DZo (DZo swapped so T1=F, T2=M)
# ===========================

rm(list=ls())
library(OpenMx)
library(dplyr)
library(readr)

# ---- Load data ----
#vcu25w_imp_06 <- read_csv("/Users/sadeghin/research-IRTAs/Amanda/from_server_data/01_tidy_data/vcu25w_imp_06.csv")
vcu25w_imp_06 <- read_csv("/Users/sadeghin/research-IRTAs/Amanda/scripts_for_paper/vcu25w_imp_paper.csv")


# ---- User inputs ----
vars <- c("rsrs_imp","rsds_imp")
nv   <- length(vars)
ntv  <- nv*2
selVars <- paste(vars, c(rep("_T1",nv), rep("_T2",nv)), sep="")

# ---- Prepare data (ensure numeric) ----
df <- vcu25w_imp_06
for(v in selVars){
  if(is.factor(df[[v]]) || is.ordered(df[[v]])){
    df[[v]] <- as.numeric(as.character(df[[v]]))
  }
}

# ---- Split groups ----
mzfData <- subset(df, (zygo=="MZ" & sex_T1==2 & sex_T2==2), select=selVars)
dzfData <- subset(df, (zygo=="DZ" & sex_T1==2 & sex_T2==2), select=selVars)
mzmData <- subset(df, (zygo=="MZ" & sex_T1==1 & sex_T2==1), select=selVars)
dzmData <- subset(df, (zygo=="DZ" & sex_T1==1 & sex_T2==1), select=selVars)

# ---- DZo: swap so T1 = female, T2 = male ----
dzoRaw <- subset(df, (zygo=="DZ" & sex_T1 != sex_T2), select=c("sex_T1","sex_T2", selVars))
dzOS <- dzoRaw
swap <- dzoRaw$sex_T1==1 & dzoRaw$sex_T2==2
if(any(swap, na.rm=TRUE)){
  tmp <- dzOS[swap, selVars[1:nv]]
  dzOS[swap, selVars[1:nv]] <- dzOS[swap, selVars[(nv+1):(2*nv)]]
  dzOS[swap, selVars[(nv+1):(2*nv)]] <- tmp
  tmpS <- dzOS$sex_T1[swap]
  dzOS$sex_T1[swap] <- dzOS$sex_T2[swap]
  dzOS$sex_T2[swap] <- tmpS
}
dzoData <- subset(dzOS, select=selVars)

cat("Group sizes (MZf, DZf, MZm, DZm, DZo):\n")
print(c(nrow(mzfData), nrow(dzfData), nrow(mzmData), nrow(dzmData), nrow(dzoData)))

# ---- Helper for lower-tri labels ----
labLowerTri <- function(prefix, n){
  labs <- matrix("", n, n)
  for(i in 1:n) for(j in 1:i) labs[i,j] <- paste0(prefix, i, j)
  labs[lower.tri(labs, diag=TRUE)]
}

# ---- Build one group model (self-contained) ----
# groupType: "MZf","DZf","MZm","DZm","DZo"
# freeR: whether to free rgA_par and rcC_par in DZo
buildGroup <- function(groupType, data, selVars, nv, freeR=FALSE){
  
  # Means: sex-specific mean vectors (shared across groups by label)
  meanF <- mxMatrix("Full", 1, nv, free=TRUE, values=0,
                    labels=paste0("meanF_", 1:nv), name="meanF")
  meanM <- mxMatrix("Full", 1, nv, free=TRUE, values=0,
                    labels=paste0("meanM_", 1:nv), name="meanM")
  
  # Cholesky paths (shared across groups by label)
  aF <- mxMatrix("Lower", nv, nv, free=TRUE, values=.3, labels=labLowerTri("aF",nv), name="aF")
  cF <- mxMatrix("Lower", nv, nv, free=TRUE, values=.2, labels=labLowerTri("cF",nv), name="cF")
  eF <- mxMatrix("Lower", nv, nv, free=TRUE, values=.4, labels=labLowerTri("eF",nv), name="eF")
  
  aM <- mxMatrix("Lower", nv, nv, free=TRUE, values=.3, labels=labLowerTri("aM",nv), name="aM")
  cM <- mxMatrix("Lower", nv, nv, free=TRUE, values=.2, labels=labLowerTri("cM",nv), name="cM")
  eM <- mxMatrix("Lower", nv, nv, free=TRUE, values=.4, labels=labLowerTri("eM",nv), name="eM")
  
  # Variance components
  Af <- mxAlgebra(aF %*% t(aF), name="Af")
  Cf <- mxAlgebra(cF %*% t(cF), name="Cf")
  Ef <- mxAlgebra(eF %*% t(eF), name="Ef")
  Pf <- mxAlgebra(Af + Cf + Ef, name="Pf")
  
  Am <- mxAlgebra(aM %*% t(aM), name="Am")
  Cm <- mxAlgebra(cM %*% t(cM), name="Cm")
  Em <- mxAlgebra(eM %*% t(eM), name="Em")
  Pm <- mxAlgebra(Am + Cm + Em, name="Pm")
  
  # Mean row and expected covariance (local to group)
  if(groupType == "MZf"){
    meanRow <- mxAlgebra(cbind(meanF, meanF), name="expMean")
    cross   <- mxAlgebra(Af + Cf, name="cross")
    expCov  <- mxAlgebra(rbind(cbind(Pf, cross),
                               cbind(t(cross), Pf)), name="expCov")
    
  } else if(groupType == "DZf"){
    meanRow <- mxAlgebra(cbind(meanF, meanF), name="expMean")
    cross   <- mxAlgebra(0.5*Af + Cf, name="cross")
    expCov  <- mxAlgebra(rbind(cbind(Pf, cross),
                               cbind(t(cross), Pf)), name="expCov")
    
  } else if(groupType == "MZm"){
    meanRow <- mxAlgebra(cbind(meanM, meanM), name="expMean")
    cross   <- mxAlgebra(Am + Cm, name="cross")
    expCov  <- mxAlgebra(rbind(cbind(Pm, cross),
                               cbind(t(cross), Pm)), name="expCov")
    
  } else if(groupType == "DZm"){
    meanRow <- mxAlgebra(cbind(meanM, meanM), name="expMean")
    cross   <- mxAlgebra(0.5*Am + Cm, name="cross")
    expCov  <- mxAlgebra(rbind(cbind(Pm, cross),
                               cbind(t(cross), Pm)), name="expCov")
    
  } else if(groupType == "DZo"){
    # DZo: T1=F, T2=M (because you swapped)
    meanRow <- mxAlgebra(cbind(meanF, meanM), name="expMean")
    
    # Scalars ONLY in DZo
    rgA <- mxMatrix("Full", 1, 1, free=freeR, values=ifelse(freeR, .8, 1),
                    labels="rgA_par", lbound=-1, ubound=1, name="rgA")
    rcC <- mxMatrix("Full", 1, 1, free=freeR, values=ifelse(freeR, .8, 1),
                    labels="rcC_par", lbound=-1, ubound=1, name="rcC")
    
    A_FM <- mxAlgebra(rgA[1,1] * (aF %*% t(aM)), name="A_FM")
    C_FM <- mxAlgebra(rcC[1,1] * (cF %*% t(cM)), name="C_FM")
    
    cross <- mxAlgebra(0.5*A_FM + C_FM, name="cross")
    expCov <- mxAlgebra(rbind(cbind(Pf, cross),
                              cbind(t(cross), Pm)), name="expCov")
    
  } else stop("Unknown groupType")
  
  dataObj <- mxData(observed=data, type="raw")
  expObj  <- mxExpectationNormal(covariance="expCov", means="expMean", dimnames=selVars)
  fitFun  <- mxFitFunctionML()
  
  if(groupType == "DZo"){
    model <- mxModel(groupType,
                     meanF, meanM,
                     aF, cF, eF, aM, cM, eM,
                     Af, Cf, Ef, Pf, Am, Cm, Em, Pm,
                     rgA, rcC, A_FM, C_FM,
                     meanRow, cross, expCov,
                     dataObj, expObj, fitFun)
  } else {
    model <- mxModel(groupType,
                     meanF, meanM,
                     aF, cF, eF, aM, cM, eM,
                     Af, Cf, Ef, Pf, Am, Cm, Em, Pm,
                     meanRow, cross, expCov,
                     dataObj, expObj, fitFun)
  }
  model
}

# ===========================
# Quantitative ACE (rgA=1, rcC=1 in DZo)
# ===========================
modelMZf <- buildGroup("MZf", mzfData, selVars, nv, freeR=FALSE)
modelDZf <- buildGroup("DZf", dzfData, selVars, nv, freeR=FALSE)
modelMZm <- buildGroup("MZm", mzmData, selVars, nv, freeR=FALSE)
modelDZm <- buildGroup("DZm", dzmData, selVars, nv, freeR=FALSE)
modelDZo <- buildGroup("DZo", dzoData, selVars, nv, freeR=FALSE)


modelACE_quant <- mxModel("ACE_quant",
                          modelMZf, modelDZf, modelMZm, modelDZm, modelDZo,
                          mxFitFunctionMultigroup(c("MZf","DZf","MZm","DZm","DZo")))

cat("\nRunning quantitative ACE (rgA=1, rcC=1)...\n")
fitACE_quant <- mxTryHard(modelACE_quant, extraTries=10, intervals=TRUE)
summary(fitACE_quant)

# ===========================
# Full sex-limitation ACE (free rgA, rcC in DZo)
# ===========================
modelMZf2 <- buildGroup("MZf", mzfData, selVars, nv, freeR=TRUE)  # freeR won't matter here
modelDZf2 <- buildGroup("DZf", dzfData, selVars, nv, freeR=TRUE)
modelMZm2 <- buildGroup("MZm", mzmData, selVars, nv, freeR=TRUE)
modelDZm2 <- buildGroup("DZm", dzmData, selVars, nv, freeR=TRUE)
modelDZo2 <- buildGroup("DZo", dzoData, selVars, nv, freeR=TRUE)

modelACE_full <- mxModel("ACE_full",
                         modelMZf2, modelDZf2, modelMZm2, modelDZm2, modelDZo2,
                         mxFitFunctionMultigroup(c("MZf","DZf","MZm","DZm","DZo")))

cat("\nRunning FULL sex-lim ACE (free rgA, rcC)...\n")
fitACE_full <- mxTryHard(modelACE_full, extraTries=10, intervals=FALSE)
summary(fitACE_full)

cat("\nCompare Quant vs Full (LRT):\n")

print(mxCompare(fitACE_full, fitACE_quant))

# Optional: show rgA, rcC estimates from full model
cat("\nFull-model rgA, rcC (DZo) estimates:\n")
print(omxGetParameters(fitACE_full)[c("rgA_par","rcC_par")])

# ============================================================
# ONE-STOP EXTRACT + PRINT EVERYTHING (A/C/E, h2/c2/e2, rA/rC/rE)
# Works with your fitted model objects: fitACE_quant / fitACE_full
# Paste this whole block after the model has successfully run.
# ============================================================
# ============================================================
# ONE-STOP EXTRACT + PRINT EVERYTHING (robust to nesting)
# Works with fitACE_quant / fitACE_full.
# ============================================================

fit <- fitACE_quant   # or fitACE_quant

# ---- recursively search model tree for an entity name and mxEvalByName it
mxEvalFind <- function(model, name) {
  # try here
  out <- try(mxEvalByName(name, model), silent=TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  # recurse into submodels
  kids <- names(model)
  for (k in kids) {
    obj <- try(model[[k]], silent=TRUE)
    if (!inherits(obj, "try-error") && inherits(obj, "MxModel")) {
      out2 <- try(mxEvalFind(obj, name), silent=TRUE)
      if (!inherits(out2, "try-error")) return(out2)
    }
  }
  stop(paste0("Could not find entity named '", name, "' anywhere in model tree."))
}

printMat <- function(x, title=NULL, digits=4) {
  if (!is.null(title)) cat("\n", title, "\n", sep="")
  print(round(x, digits))
}

# ---- pull A/C/E/P (female + male) wherever they live
Af <- mxEvalFind(fit, "Af")
Cf <- mxEvalFind(fit, "Cf")
Ef <- mxEvalFind(fit, "Ef")
Pf <- mxEvalFind(fit, "Pf")

Am <- mxEvalFind(fit, "Am")
Cm <- mxEvalFind(fit, "Cm")
Em <- mxEvalFind(fit, "Em")
Pm <- mxEvalFind(fit, "Pm")

# ---- standardized components
h2_f <- Af / Pf; c2_f <- Cf / Pf; e2_f <- Ef / Pf
h2_m <- Am / Pm; c2_m <- Cm / Pm; e2_m <- Em / Pm

# ---- correlations
rA_f <- cov2cor(Af); rC_f <- cov2cor(Cf); rE_f <- cov2cor(Ef)
rA_m <- cov2cor(Am); rC_m <- cov2cor(Cm); rE_m <- cov2cor(Em)

# ---- cross-sex scalars (these may only exist in DZo submodel, depending on how you built it)
rgA_val <- try(mxEvalFind(fit, "rgA_unique_4n1"), silent=TRUE)
rcC_val <- try(mxEvalFind(fit, "rcC_unique_4n1"), silent=TRUE)

# ---- means
meanF_val <- mxEvalFind(fit, "meanF")
meanM_val <- mxEvalFind(fit, "meanM")

cat("\n==================== MODEL SUMMARY ====================\n")
cat("Model name:", fit$name, "\n")
cat("minus2LL:", fit$output$Minus2LogLikelihood, " df:", fit$output$degreesOfFreedom, "\n")
cat("AIC:", fit$output$AIC, "\n")

cat("\n==================== FEMALE ====================\n")
printMat(Af, "Af (female)")
printMat(Cf, "Cf (female)")
printMat(Ef, "Ef (female)")
printMat(Pf, "Pf (female)")
printMat(h2_f, "h2_f = Af/Pf (female)")
printMat(c2_f, "c2_f = Cf/Pf (female)")
printMat(e2_f, "e2_f = Ef/Pf (female)")
printMat(rA_f, "rA_f = cor(Af) (female)")
printMat(rC_f, "rC_f = cor(Cf) (female)")
printMat(rE_f, "rE_f = cor(Ef) (female)")

cat("\n==================== MALE ====================\n")
printMat(Am, "Am (male)")
printMat(Cm, "Cm (male)")
printMat(Em, "Em (male)")
printMat(Pm, "Pm (male)")
printMat(h2_m, "h2_m = Am/Pm (male)")
printMat(c2_m, "c2_m = Cm/Pm (male)")
printMat(e2_m, "e2_m = Em/Pm (male)")
printMat(rA_m, "rA_m = cor(Am) (male)")
printMat(rC_m, "rC_m = cor(Cm) (male)")
printMat(rE_m, "rE_m = cor(Em) (male)")

cat("\n==================== CROSS-SEX SCALARS ====================\n")
if (!inherits(rgA_val, "try-error")) cat("rgA_unique_4n1 =", as.numeric(rgA_val), "\n") else cat("rgA_unique_4n1 not found (ok if not in model)\n")
if (!inherits(rcC_val, "try-error")) cat("rcC_unique_4n1 =", as.numeric(rcC_val), "\n") else cat("rcC_unique_4n1 not found (ok if not in model)\n")

cat("\n==================== MEANS ====================\n")
printMat(meanF_val, "meanF (female means)")
printMat(meanM_val, "meanM (male means)")

cat("\nDone.\n")



# ============================================================
# AE model (drop C from quantitative ACE)
# ============================================================

# Start from the quantitative ACE model
modelAE_quant <- mxModel(fitACE_quant, name="AE_quant")

# Get parameter names
pars <- names(omxGetParameters(modelAE_quant))

# Identify all C-path parameters (cF and cM)
c_pars <- pars[grep("^cF|^cM", pars)]

cat("Fixing these C parameters to zero:\n")
print(c_pars)

# Fix them to zero
for(p in c_pars){
  modelAE_quant <- omxSetParameters(modelAE_quant,
                                    labels=p,
                                    free=FALSE,
                                    values=0)
}

# Run AE model
cat("\nRunning AE model (C fixed to zero)...\n")
fitAE_quant <- mxTryHard(modelAE_quant, extraTries=10, intervals=TRUE)
summary(fitAE_quant)

fitAE_quant$output$confidenceIntervals

# Compare AE vs ACE_quant
cat("\nCompare AE vs ACE_quant:\n")
print(mxCompare(fitACE_quant, fitAE_quant))



#**********
#*# ============================================================
# Quantitative CE model: drop A from the quantitative ACE model
# ============================================================

# Start from the quantitative ACE model
modelCE_quant <- mxModel(fitACE_quant, name = "CE_quant")

# Get parameter names
pars <- names(omxGetParameters(modelCE_quant))

# Identify all A-path parameters (aF and aM)
a_pars <- pars[grep("^(aF|aM)", pars)]

cat("Fixing these A parameters to zero:\n")
print(a_pars)

# Fix them to zero
for (p in a_pars) {
  modelCE_quant <- omxSetParameters(
    modelCE_quant,
    labels = p,
    free   = FALSE,
    values = 0
  )
}

# Run CE model
cat("\nRunning CE model (A fixed to zero)...\n")
fitCE_quant <- mxTryHard(modelCE_quant, extraTries = 10, intervals = FALSE)

summary(fitCE_quant)

# Compare CE vs ACE_quant
cat("\nCompare CE vs ACE_quant:\n")
print(mxCompare(fitACE_quant, fitCE_quant))

# ============================================================
# E-only model: drop A and C from the quantitative ACE model
# ============================================================

# Start from the quantitative ACE model
modelE_quant <- mxModel(fitACE_quant, name = "E_quant")

# Get all parameter names
pars <- names(omxGetParameters(modelE_quant))

# Identify A and C path parameters
ac_pars <- pars[grep("^(aF|aM|cF|cM)", pars)]

cat("Fixing these A and C parameters to zero:\n")
print(ac_pars)

# Fix them to zero
for (p in ac_pars) {
  modelE_quant <- omxSetParameters(
    modelE_quant,
    labels = p,
    free   = FALSE,
    values = 0
  )
}

# Run E model
cat("\nRunning E-only model (A and C fixed to zero)...\n")
fitE_quant <- mxTryHard(modelE_quant, extraTries = 10, intervals = FALSE)

summary(fitE_quant)

# Compare E vs ACE_quant
cat("\nCompare E vs ACE_quant:\n")
print(mxCompare(fitACE_quant, fitE_quant))

#*************

#***********
#*
#*
#
fit <- fitAE_quant   # or fitACE_quant

# ---- recursively search model tree for an entity name and mxEvalByName it
mxEvalFind <- function(model, name) {
  # try here
  out <- try(mxEvalByName(name, model), silent=TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  # recurse into submodels
  kids <- names(model)
  for (k in kids) {
    obj <- try(model[[k]], silent=TRUE)
    if (!inherits(obj, "try-error") && inherits(obj, "MxModel")) {
      out2 <- try(mxEvalFind(obj, name), silent=TRUE)
      if (!inherits(out2, "try-error")) return(out2)
    }
  }
  stop(paste0("Could not find entity named '", name, "' anywhere in model tree."))
}

printMat <- function(x, title=NULL, digits=4) {
  if (!is.null(title)) cat("\n", title, "\n", sep="")
  print(round(x, digits))
}

# ---- pull A/C/E/P (female + male) wherever they live
Af <- mxEvalFind(fit, "Af")
Cf <- mxEvalFind(fit, "Cf")
Ef <- mxEvalFind(fit, "Ef")
Pf <- mxEvalFind(fit, "Pf")

Am <- mxEvalFind(fit, "Am")
Cm <- mxEvalFind(fit, "Cm")
Em <- mxEvalFind(fit, "Em")
Pm <- mxEvalFind(fit, "Pm")

# ---- standardized components
h2_f <- Af / Pf; c2_f <- Cf / Pf; e2_f <- Ef / Pf
h2_m <- Am / Pm; c2_m <- Cm / Pm; e2_m <- Em / Pm

# ---- correlations
rA_f <- cov2cor(Af); rC_f <- cov2cor(Cf); rE_f <- cov2cor(Ef)
rA_m <- cov2cor(Am); rC_m <- cov2cor(Cm); rE_m <- cov2cor(Em)

# ---- cross-sex scalars (these may only exist in DZo submodel, depending on how you built it)
rgA_val <- try(mxEvalFind(fit, "rgA_unique_4n1"), silent=TRUE)
rcC_val <- try(mxEvalFind(fit, "rcC_unique_4n1"), silent=TRUE)

# ---- means
meanF_val <- mxEvalFind(fit, "meanF")
meanM_val <- mxEvalFind(fit, "meanM")

cat("\n==================== MODEL SUMMARY ====================\n")
cat("Model name:", fit$name, "\n")
cat("minus2LL:", fit$output$Minus2LogLikelihood, " df:", fit$output$degreesOfFreedom, "\n")
cat("AIC:", fit$output$AIC, "\n")

cat("\n==================== FEMALE ====================\n")
printMat(Af, "Af (female)")
printMat(Cf, "Cf (female)")
printMat(Ef, "Ef (female)")
printMat(Pf, "Pf (female)")
printMat(h2_f, "h2_f = Af/Pf (female)")
printMat(c2_f, "c2_f = Cf/Pf (female)")
printMat(e2_f, "e2_f = Ef/Pf (female)")
printMat(rA_f, "rA_f = cor(Af) (female)")
printMat(rC_f, "rC_f = cor(Cf) (female)")
printMat(rE_f, "rE_f = cor(Ef) (female)")

cat("\n==================== MALE ====================\n")
printMat(Am, "Am (male)")
printMat(Cm, "Cm (male)")
printMat(Em, "Em (male)")
printMat(Pm, "Pm (male)")
printMat(h2_m, "h2_m = Am/Pm (male)")
printMat(c2_m, "c2_m = Cm/Pm (male)")
printMat(e2_m, "e2_m = Em/Pm (male)")
printMat(rA_m, "rA_m = cor(Am) (male)")
printMat(rC_m, "rC_m = cor(Cm) (male)")
printMat(rE_m, "rE_m = cor(Em) (male)")

cat("\n==================== CROSS-SEX SCALARS ====================\n")
if (!inherits(rgA_val, "try-error")) cat("rgA_unique_4n1 =", as.numeric(rgA_val), "\n") else cat("rgA_unique_4n1 not found (ok if not in model)\n")
if (!inherits(rcC_val, "try-error")) cat("rcC_unique_4n1 =", as.numeric(rcC_val), "\n") else cat("rcC_unique_4n1 not found (ok if not in model)\n")

cat("\n==================== MEANS ====================\n")
printMat(meanF_val, "meanF (female means)")
printMat(meanM_val, "meanM (male means)")

cat("\nDone.\n")


############
# ---------- helper to safely mxEval at the right level ----------
mxEvalSmart <- function(expr, fit) {
  # try top model first
  out <- try(mxEval(expr, fit), silent=TRUE)
  if (!inherits(out, "try-error")) return(out)
  
  # try common submodels
  for (nm in c("MZf","DZf","MZm","DZm","DZo")) {
    if (!is.null(fit[[nm]])) {
      out <- try(mxEval(expr, fit[[nm]]), silent=TRUE)
      if (!inherits(out, "try-error")) return(out)
    }
  }
  stop("Couldn't mxEval expression at any model level: ", deparse(substitute(expr)))
}

# ---------- decompose covariance for one sex ----------
cov_decomp <- function(fit, sex=c("f","m")) {
  sex <- match.arg(sex)
  
  if (sex=="f") {
    A <- mxEvalSmart(Af, fit)
    C <- mxEvalSmart(Cf, fit)
    E <- mxEvalSmart(Ef, fit)
    P <- mxEvalSmart(Pf, fit)
    tag <- "Female"
  } else {
    A <- mxEvalSmart(Am, fit)
    C <- mxEvalSmart(Cm, fit)
    E <- mxEvalSmart(Em, fit)
    P <- mxEvalSmart(Pm, fit)
    tag <- "Male"
  }
  
  covA12 <- A[1,2]; covC12 <- C[1,2]; covE12 <- E[1,2]; covP12 <- P[1,2]
  sd1 <- sqrt(P[1,1]); sd2 <- sqrt(P[2,2])
  
  out <- data.frame(
    Sex = tag,
    CovP12 = covP12,
    CovA12 = covA12,
    CovC12 = covC12,
    CovE12 = covE12,
    PropA_of_CovP12 = covA12 / covP12,
    PropC_of_CovP12 = covC12 / covP12,
    PropE_of_CovP12 = covE12 / covP12,
    rP12 = covP12 / (sd1*sd2),
    rA12 = covA12 / (sd1*sd2),
    rC12 = covC12 / (sd1*sd2),
    rE12 = covE12 / (sd1*sd2)
  )
  return(out)
}

# ---------- run it ----------
# set fit to whichever model you want to report (ACE_full, ACE_quant, etc.)
fit <- fitAE_quant   # <- change if needed

print(rbind(
  cov_decomp(fit, "f"),
  cov_decomp(fit, "m")
))


############
### Added May 8
#######
# ============================================================
# Cross-sex comparison set:
# quant, rgA-only, rcC-only, full
# ============================================================


# =============================
# rgA-only model
# =============================
model_rgAonly <- mxModel(fitACE_quant, name="ACE_rgAonly")

model_rgAonly <- omxSetParameters(model_rgAonly,
                                  labels="rgA_par",
                                  free=TRUE, values=0.8,
                                  lbound=-1, ubound=1)

model_rgAonly <- omxSetParameters(model_rgAonly,
                                  labels="rcC_par",
                                  free=FALSE, values=1)

fit_rgAonly <- mxTryHard(model_rgAonly, extraTries=10)

# =============================
# rcC-only model
# =============================
model_rcConly <- mxModel(fitACE_quant, name="ACE_rcConly")

model_rcConly <- omxSetParameters(model_rcConly,
                                  labels="rcC_par",
                                  free=TRUE, values=0.8,
                                  lbound=-1, ubound=1)

model_rcConly <- omxSetParameters(model_rcConly,
                                  labels="rgA_par",
                                  free=FALSE, values=1)

fit_rcConly <- mxTryHard(model_rcConly, extraTries=10)

fit_row <- function(fit, model_name) {
  if (inherits(fit, "try-error")) {
    return(data.frame(
      model = model_name,
      ep = NA,
      minus2LL = NA,
      df = NA,
      AIC = NA
    ))
  }
  
  s <- summary(fit)
  
  data.frame(
    model = model_name,
    ep = length(omxGetParameters(fit)),
    minus2LL = s$Minus2LogLikelihood,
    df = s$degreesOfFreedom,
    AIC = s$AIC
  )
}

# Model fit table
fit_table <- rbind(
  fit_row(fitACE_quant, "ACE_quant"),
  fit_row(fitAE_quant,  "AE_quant"),
  fit_row(fitCE_quant,  "CE_quant"),
  fit_row(fitE_quant,  "E_quant"),
  fit_row(fit_rgAonly,  "rgA_only"),
  fit_row(fit_rcConly,  "rcC_only"),
  fit_row(fitACE_full,  "ACE_full")
)

print(fit_table)

# Nested likelihood-ratio tests
lrt_AE_vs_ACE   <- mxCompare(fitACE_quant, fitAE_quant)
lrt_rgA_vs_ACE  <- mxCompare(fit_rgAonly,fitACE_quant)
lrt_rcC_vs_ACE  <- mxCompare(fit_rcConly,fitACE_quant)
lrt_full_vs_ACE <- mxCompare(fitACE_full,fitACE_quant)

cat("\n--- LRT: AE vs ACE_quant ---\n")
print(lrt_AE_vs_ACE)

cat("\n--- LRT: rgA_only vs ACE_quant ---\n")
print(lrt_rgA_vs_ACE)

cat("\n--- LRT: rcC_only vs ACE_quant ---\n")
print(lrt_rcC_vs_ACE)

cat("\n--- LRT: full vs ACE_quant ---\n")
print(lrt_full_vs_ACE)

# Cross-sex scalar estimates
cat("\n--- Cross-sex scalar estimates ---\n")
cat("rgA_only model:\n")
print(omxGetParameters(fit_rgAonly)[c("rgA_par")])

cat("\nrcC_only model:\n")
print(omxGetParameters(fit_rcConly)[c("rcC_par")])

cat("\nfull model:\n")
print(omxGetParameters(fitACE_full)[c("rgA_par", "rcC_par")])

mxCompare( vcu_fitSAT, nested <- list( fitACE_quant, fitAE_quant, fitCE_quant, fitE_quant) )
mxCompare( vcu_fitSAT, nested <- list( fitACE_full, fitAE, fitCE, fitE) )






