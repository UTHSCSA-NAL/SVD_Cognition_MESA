# Load required libraries
libraries <- c(
  'psych',
  'GPArotation',
  'lavaan',
  'semTools',
  'semPlot',
  'semptools',
  'data.table')

for (p in libraries) {
  if (!requireNamespace(p)) {
    install.packages(p)
  }}

lapply(libraries, require, character.only = TRUE)

rm(list = c("libraries", "p"))

# Name dictionary:

# mri6 = MRI dataset MESA exam 6
# mrivar = MRI marker variable vector
# cogvar = cognitive score variable vector
# dtabs = time interval between MRI scan and cognitive assessment
# dtord = order of MRI scan and cognitive assessment completion
# z_bgpvsc = basal ganglia PVS count
# z_thpvsc = thalamus PVS count
# z_pvwml = periventricular WMH volume
# z_otwml = subcortical WMH volume
# z_FA_MUSE_604 = mean white matter fractional anisotropy
# z_TR_MUSE_604 = mean white matter trace
# z_MUSE_601 = total gray matter volume
# z_MUSE_702 = total intracranial volume
# GLOBAL_COG = global cognitive score
# frci086c = Framingham risk score for global cardiovascular disease (FRS) at MESA exam 6
# site = dummy variable vector for study sites with "University of Minnesota" as reference
# race = dummy variable vector for race/ethnicity with "White" as reference
# edu = dummy variable vector for education with "Highschool or lower" as reference
# lng = dummy variable vector for language of cognitive testing with "English" as reference

#//---------------------------------------------------------------------------START OF SOURCE CODE---------------------------------------------------------------------------//

#//--------------------------------------------ANALYSES 2.5.1: MRI MARKER DIMENSIONS AND ASSOCIATIONS WITH COGNITIVE PERFORMANCE--------------------------------------------//

# Create a mixed correlation matrix with Pearson (between continuous), tetrachoric (between dichotomous), and biserial (between continuous and dichotomous) correlations
  # cmb = prefix of microbleed variables
mix_cor <- psych::mixedCor(data = mri6[,mrivar],
                    c = grep("cmb", colnames(mri6[,mrivar]), invert = TRUE),
                    d = grep("cmb", colnames(mri6[,mrivar])), method = "pearson",
                    smooth = FALSE, correct = FALSE, global = FALSE)

# Correlation matrix diagnostics
mix_cor$rho # Correlation matrix
psych::smc(mix_cor$rho) # Calculate squared multiple correlations to assess for multicollinearity
det(mix_cor$rho) # Calculate matrix determinant to assess for singularity of correlation matrix
psych::KMO(mix_cor$rho) # Calculate KMO test

# Determine number of factors to be extracted based on the Very Simple Structure and Minimum Average Partial criteria
vss <- psych::VSS(x = mix_cor$rho, n = 8, rotate = "oblimin", fm = "minres",
                  n.obs = dim(mri6)[1])
vss

# Plot MAP criterion
ggplot(mapping = aes(y = vss$map, x = c(1:length(vss$map)))) + theme_linedraw() +
  geom_line(linewidth = 0.6) + ggtitle("Minimum Average Partial") +
  ylab("Minimum Average Partial") + xlab("Number of Components") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = 7),
        axis.title.x = element_text(size = 12, vjust = -5),
        axis.title.y = element_text(size = 12, vjust = 8)) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(round(min(vss$map),2), round(max(vss$map),2), 0.01)) +
  scale_x_continuous(breaks = seq(1, length(vss$map), 1)) +
  theme(plot.margin = margin(1.5,1,1.5,1, "cm"))

# Create a function that runs an EFA on the correlation matrix iteratively, removing the item with the lowest communality from the matrix at each iteration, until all remaining items have communalities >=0.2
# The function saves the final solution and also calculates and saves factor scores
# We also include the option for extension analysis
efa_fun <- function(var, data, factors, rotation = "oblimin", fm = "minres", scores = "regression",
                    com = 0.2, ext = FALSE, ext.proc = "D"){
  var_red <- var
  repeat{
    if(any(grepl("cmb", var_red))){
      mat <- psych::mixedCor(data = data[,var_red],
                      c = grep("cmb", var_red, value = TRUE, invert = TRUE),
                      d = grep("cmb", var_red, value = TRUE),
                      smooth = FALSE, correct = FALSE, global = FALSE)$rho
    }
    else{
      mat <- mat <- cor(x = data[,var_red], method = "pearson")}
    fac <- psych::fa(r = mat, nfactors = factors, n.obs = dim(data)[1],
                     rotate = rotation, fm = fm)
    print(round(sort(fac$communality, decreasing = FALSE), 2))
    if(all(fac$communality >= com)) break
    var_red <- var_red[!var_red %in% names(sort(fac$communality, decreasing = FALSE))[1]]
  }
  fac_scores <- as.data.frame(psych::factor.scores(x = data[,var_red],
                                                   f = fac, method = scores)$scores)
  if(ext == TRUE){
    if(ext.proc == "D"){proc <- TRUE}
    if(ext.proc == "G"){proc <- FALSE}
    
    mat_ext <- psych::mixedCor(data = data[,var],
                        c = grep("cmb", var, value = TRUE, invert = TRUE),
                        d = grep("cmb", var, value = TRUE),
                        smooth = FALSE, correct = FALSE, global = FALSE)$rho
    
    fac_ext <- psych::fa.extend(r = mat_ext, nfactors = factors,
                         ov = which(colnames(mat_ext) %in% var_red),
                         ev = which(!colnames(mat_ext) %in% var_red), n.obs = dim(data)[1],
                         correct = proc, rotate = rotation, fm = fm)
    
    assign(paste("fac_ext", deparse(substitute(var)), sep = "_"), fac_ext, envir = globalenv())
  }
  assign(paste("fac", deparse(substitute(var)), sep = "_"), fac, envir = globalenv())
  assign(paste("scor", deparse(substitute(var)), sep = "_"), fac_scores, envir = globalenv())
}

# Run function for all SVD marker variables with oblimin rotation
efa_fun(var = mrivar, data = mri6, factors = 2, rotation = "oblimin",
        fm = "minres", scores = "regression", com = 0.2)

# Print EFA solution with oblimin rotation
print(mget(ls(pattern = "fac_mrivar")))

fac_mrivar$Phi # interfactor correlation matrix

# Run function for all SVD marker variables with varimax rotation 
# Also perform extension analysis using the Dwyer's factor extension method
efa_fun(var = mrivar, data = mri6, factors = 2, rotation = "varimax",
        fm = "minres", ext = TRUE, scores = "regression",
        ext.proc = "D")

# Print EFA solution with varimax rotation
print(mget(ls(pattern = "fac_mrivar")))

# Print results of extension analysis
print(mget(ls(pattern = "fac_ext")))

# Merge computed factor scores with main dataset
mri6$MR1_r <- scor_mrivar$MR1
mri6$MR2_r <- scor_mrivar$MR2

# Associations of the factor scores with cognitive performance
# Create a loop GLM function with 95% profile-likelihood-based CIs, and FDR control for multiple predictors (MRI markers)
glm_mod <- function(outc, pred, adj = NULL, data, p.adj = FALSE, c.int = FALSE) {
  for (i in seq_along(pred)){
    fml <- paste0(outc, "~", paste(adj, collapse = "+"), "+", pred[i])
    fml <- as.formula(fml)
    fit <- glm(formula = fml, data = data, family = gaussian)
    cat("Predictor:", pred[i], "; Participants included:", nobs(fit), "\n")
    mat <- as.data.frame(summary.glm(fit)$coefficients[,c("Estimate", "Std. Error", "Pr(>|t|)")])
    mat <- mat[pred[i], , drop = FALSE]
    if(c.int == TRUE){
      mat <- cbind(mat, t(confint(profile(fit), parm = pred[i])))}
    assign(pred[i], mat)
    if(all(pred %in% ls(pattern = paste(pred, collapse = "|")))){
      
      predres <- do.call(bind_rows, mget(pred))
      rm(list = pred)
      if(p.adj == TRUE){
        predres$Adj.p <- p.adjust(predres[,"Pr(>|t|)"], method = "BH")
        predres$Sig <- ifelse(predres$Adj.p <= 0.05, "YES", "")
        predres <- predres[,!colnames(predres) %in% "Pr(>|t|)"]}
      if(p.adj == FALSE){
        predres$Sig <- ifelse(predres[,"Pr(>|t|)"] <= 0.05, "YES", "")}
      predres <- predres %>% mutate(across(where(is.numeric), function(x) {
        case_when(abs(x) >= 0.001 ~ formatC(x, format = "f", digits = 3), TRUE ~ formatC(x, format = "e", digits = 2))}))
      if(c.int == TRUE){
        predres$ci <- paste(predres$`2.5 %`, predres$`97.5 %`, sep = ", ")
      }
      colnames(predres) <- paste(outc, colnames(predres), sep = "_")
      predres <- as.data.frame(t(predres))
      assign(paste0("glm_mod_", outc), predres, envir = globalenv())}
  }}

# Create loop function that runs the GLM function across all cognitive scores
glm_mod_loop <- function(out, pred, adj = NULL, data, p.adj = FALSE, c.int = FALSE){
  for (j in seq_along(out))
  {
    cat("\n> Analysis ", j, " of ", length(out), " <", "\nOutcome: ", out[j], "\n")
    glm_mod(outc = paste0(out[j]), pred = pred, adj = adj, data = data, p.adj = p.adj, c.int = c.int)
    if(all(paste0("glm_mod_", out) %in% ls(pattern = paste(out, collapse = "|"), envir = globalenv()))){
      glmres <- do.call(bind_rows, mget(paste0("glm_mod_", out), envir = globalenv()))
      assign("glmres", glmres, envir = globalenv())
      rm(list = ls(pattern = paste(paste0("glm_mod_", out), collapse = "|"), envir = globalenv()), envir = globalenv())
    }}}

glm_mod_loop(out = cogvar, pred = c("MR1_r", "MR2_r"), p.adj = TRUE, c.int = TRUE, data = mri6,
             adj = c("AGE", "SEX_M", "z_MUSE_702", "dtabs", "dtord", race, site, edu, lng))

#//--------------------------------------------ANALYSES 2.5.2: STRUCTURAL EQUATION MODELING--------------------------------------------//
# Measurement model of SVD
model <- paste0(
  '# latent variable, measurement model
   svd =~ NA*z_bgpvsc + z_thpvsc + z_pvwml + z_otwml + z_FA_MUSE_604 + z_TR_MUSE_604
   
   # variances
   svd ~~ 1*svd
   
   # residual covariances
   z_bgpvsc ~~ z_thpvsc
   z_pvwml ~~ z_otwml
   z_FA_MUSE_604 ~~ z_TR_MUSE_604')

# Compute CFA model with mean-and-variance corrected chi-square (simple second-order correction) and robust standard errors
fit <- lavaan::sem(model, data = mri6, estimator = "ML", test = "scaled.shifted", se = "robust.sem")
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Model inspection
lavInspect(fit, "theta") # residual covariance matrix
lavInspect(fit, "sampstat") # sample covariance matrix
lavInspect(fit, "implied") # model-implied covariance matrix
lavInspect(fit, "resid")$cov # difference between observed and model-implied covariance matrix (unstandardized and unscaled model residuals)
lavResiduals(fit, type = "cor.bentler")$cov # unstandardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
lavResiduals(fit, type = "cor.bentler")$cov.z # standardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
modindices(fit, sort. = TRUE, standardized = TRUE) # sorted (from largest to smallest) model modification indices

# SEM model and Path diagram with SVD latent variable and gray matter volume as mediators in the relationship of age with global cognition
model <- paste0(
  '# latent variable, measurement model
   svd =~ NA*z_bgpvsc + z_thpvsc + z_pvwml + z_otwml + z_FA_MUSE_604 + z_TR_MUSE_604
   
   # variances
   svd ~~ 1*svd
   
   # residual covariances
   z_bgpvsc ~~ z_thpvsc
   z_pvwml ~~ z_otwml
   z_FA_MUSE_604 ~~ z_TR_MUSE_604
   
   # regressions, structural model
   svd ~ a * AGE + SEX_M +', paste(c(site, race, edu), collapse = "+"), '+ z_MUSE_702
   z_MUSE_601 ~ b * AGE + c * svd + SEX_M +', paste(c(site, race, edu), collapse = "+"), '+ z_MUSE_702
   GLOBAL_COG ~ d * svd + e * z_MUSE_601 + f * AGE + SEX_M +', paste(c(site, race, edu), collapse = "+"), '+ z_MUSE_702
   
   # direct, indirect and total effects
   # direct effect
   direct := f
   # indirect effect svd
   in_svd := a * d
   in_svd_gm := a * c * e
   indirect_svd := in_svd + in_svd_gm
   # indirect effect GM
   indirect_gm := b * e
   # total indirect
   indirect_total := indirect_svd + indirect_gm
   # total effect
   total := direct + indirect_total')

# Compute SEM model with mean-and-variance corrected chi-square (simple second-order correction) and robust standard errors
fit <- lavaan::sem(model, data = mri6, estimator = "ML", test = "scaled.shifted", se = "robust.sem")
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Model inspection
lavInspect(fit, "theta") # residual covariance matrix
lavInspect(fit, "sampstat") # sample covariance matrix
lavInspect(fit, "implied") # model-implied covariance matrix
lavInspect(fit, "resid")$cov # difference between observed and model-implied covariance matrix (unstandardized and unscaled model residuals)
lavResiduals(fit, type = "cor.bentler")$cov # unstandardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
lavResiduals(fit, type = "cor.bentler")$cov.z # standardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
modindices(fit, sort. = TRUE, standardized = TRUE) # sorted (from largest to smallest) model modification indices

# Bootstrap solution with 5000 samples
fit <- lavaan::sem(model, data = mri6, std.lv = TRUE, estimator = "ML", test = "scaled.shifted",
                   se = "bootstrap", bootstrap = 5000, parallel = "snow",
                   ncpus = (parallel::detectCores() -1), iseed = 2023)

# Bootstrap bias-corrected 95% CIs for direct, indirect, and total effects
boot_ci <- lavaan::parameterEstimates(fit, ci = TRUE, boot.ci.type = "bca.simple",
                                      level = 0.95, standardized = TRUE)
boot_ci[boot_ci$op == ":=", c("lhs", "est", "pvalue", "ci.lower", "ci.upper")]

# Path diagram
laymat <- semptools::layout_matrix(z_b =    c(1, 2),
                                   z_th =   c(1, 6),
                                   z_p =    c(1, 10),
                                   z_tw =   c(1, 14),
                                   z_F =    c(1, 18),
                                   z_T =    c(1, 22),
                                   z_M =    c(8, 12),
                                   svd =    c(4, 12),
                                   GLO =    c(6, 23),
                                   AGE =    c(6, 1))

plt <- semptools::drop_nodes(
  object = semPlot::semPlotModel(fit),
  nodes = c(race, site, edu, "SEX_M", "z_MUSE_702"))

plt <- semPlot::semPaths(plt, whatLabels = "std", layout = laymat, residuals = FALSE, edge.label.cex = 1.2)

plt <- semptools::set_edge_label_position(plt, c("z_MUSE_601 ~ svd" = 0.2))

plt <- semptools::change_node_label(plt,
                         c(svd = "SVD",
                           z_b = "Basal\nganglia\nPVS",
                           z_th = "Thalamus\nPVS",
                           z_p = "Perive-\nntricular\nWMH",
                           z_tw = "Subco-\nrtical\nWMH",
                           z_F = "WM\nFractional\nAnisotropy",
                           z_T = "WM\nTrace",
                           z_M = "GM\nVolume",
                           AGE = "Age",
                           GLO = "Global\nCognition"),
                         label.cex = 2)

plt <- semptools::mark_sig(semPaths_plot = plt, object = fit, alphas = c('--' = 0.05, '*' = 0.01, '**' = 0.001))

plt$graphAttributes$Edges$curve <- c(0, 0, 0, 0, 0, 0, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0)

plt$graphAttributes$Edges$label.margin <- c(0, 0, 0, 0, 0, 0, -0.02, -0.02, -0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0)

plot(plt)

# SEM model and Path diagram with SVD latent variable and gray matter volume as mediators in the relationship of FRS with global cognition
model <- paste0(
  '# latent variable, measurement model
   svd =~ NA*z_bgpvsc + z_thpvsc + z_pvwml + z_otwml + z_FA_MUSE_604 + z_TR_MUSE_604
   
   # variances
   svd ~~ 1*svd
   
   # residual covariances
   z_bgpvsc ~~ z_thpvsc
   z_pvwml ~~ z_otwml
   z_FA_MUSE_604 ~~ z_TR_MUSE_604

   # regressions, structural model
   svd ~ a * frci086c +', paste(c(site, race, edu), collapse = "+"), '+ z_MUSE_702
   z_MUSE_601 ~ b * frci086c + c * svd +', paste(c(site, race, edu), collapse = "+"), '+ z_MUSE_702
   GLOBAL_COG ~ d * svd + e * z_MUSE_601 + f * frci086c +', paste(c(site, race, edu), collapse = "+"), '+ z_MUSE_702

   # direct, indirect and total effects
   direct := f
   # indirect effect svd
   in_svd := a * d
   in_svd_gm := a * c * e
   indirect_svd := in_svd + in_svd_gm
   # indirect effect GM
   indirect_gm := b * e
   # total indirect
   indirect_total := indirect_svd + indirect_gm
   # total effect
   total := direct + indirect_total')

# Compute SEM model with mean-and-variance corrected chi-square (simple second-order correction) and robust standard errors
fit <- lavaan::sem(model, data = mri6, estimator = "ML", test = "scaled.shifted", se = "robust.sem")
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Model inspection
lavInspect(fit, "theta") # residual covariance matrix
lavInspect(fit, "sampstat") # sample covariance matrix
lavInspect(fit, "implied") # model-implied covariance matrix
lavInspect(fit, "resid")$cov # difference between observed and model-implied covariance matrix (unstandardized and unscaled model residuals)
lavResiduals(fit, type = "cor.bentler")$cov # unstandardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
lavResiduals(fit, type = "cor.bentler")$cov.z # standardized model residuals after transformation to correlation matrix and rescaling (by dividing the elements by the square roots of the corresponding variances of the observed covariance matrix)
modindices(fit, sort. = TRUE, standardized = TRUE) # sorted (from largest to smallest) model modification indices

# Bootstrap solution with 5000 samples
fit <- lavaan::sem(model, data = mri6, std.lv = TRUE, estimator = "ML", test = "scaled.shifted",
                   se = "bootstrap", bootstrap = 5000, parallel = "snow",
                   ncpus = (parallel::detectCores() -1), iseed = 2023)

# Bootstrap bias-corrected 95% CIs for direct, indirect, and total effects
boot_ci <- lavaan::parameterEstimates(fit, ci = TRUE, boot.ci.type = "bca.simple",
                                      level = 0.95, standardized = TRUE)
boot_ci[boot_ci$op == ":=", c("lhs", "est", "pvalue", "ci.lower", "ci.upper")]

# Path diagram
laymat <- semptools::layout_matrix(z_b =    c(1, 2),
                                   z_th =   c(1, 6),
                                   z_p =    c(1, 10),
                                   z_tw =   c(1, 14),
                                   z_F =    c(1, 18),
                                   z_T =    c(1, 22),
                                   z_M =    c(8, 12),
                                   svd =    c(4, 12),
                                   GLO =    c(6, 23),
                                   f08 =    c(6, 1))

plt <- semptools::drop_nodes(
  object = semPlot::semPlotModel(fit),
  nodes = c(race, site, edu, "z_MUSE_702"))

plt <- semPlot::semPaths(plt, layout = laymat, whatLabels = "std", residuals = FALSE, edge.label.cex = 1.2)

plt <- semptools::set_edge_label_position(plt, c("z_MUSE_601 ~ svd" = 0.2))

plt <- semptools::change_node_label(plt,
                         c(svd = "SVD",
                           z_b = "Basal\nganglia\nPVS",
                           z_th = "Thalamus\nPVS",
                           z_p = "Perive-\nntricular\nWMH",
                           z_tw = "Subco-\nrtical\nWMH",
                           z_F = "WM\nFractional\nAnisotropy",
                           z_T = "WM\nTrace",
                           z_M = "GM\nVolume",
                           f08 = "FRS",
                           GLO = "Global\nCognition"),
                         label.cex = 2)

plt <- semptools::mark_sig(semPaths_plot = plt, object = fit, alphas = c('--' = 0.05, '*' = 0.01, '**' = 0.001))

plt$graphAttributes$Edges$curve <- c(0, 0, 0, 0, 0, 0, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0)

plt$graphAttributes$Edges$label.margin <- c(0, 0, 0, 0, 0, 0, -0.02, -0.02, -0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0)

plot(plt)

# Create a function that recomputes the structural part of SEM while substituting the latent SVD variable with each of its individual indicators
  # We include options for both bootstrap as well as Monte Carlo-derived confidence intervals (for faster computation)
  # We include option for inclusion of R-Square of the dependent variable (i.e., global cognitive score)
sem_mod <- function(ind, cor = NULL, out, pred, adj, data, c.int = "montCI", rep = 5000, perc = FALSE, r.sq = FALSE, sig = FALSE) {
  
  str_dec <- function(x) {y <- stringr::str_split_1(x, pattern = "(\\W+\\W+\\W+|\\W+\\W+|\\W+)")
  z <- stringr::str_subset(y, pattern = "(\\w+)")
  return(z)}
  
  if(length(str_dec(ind)) > 1){
    base_model <- paste0('lv =~', paste(ind), '\n',
                         paste(cor, collapse = "\n"), '\n',
                         'lv ~ a *', pred, '+', paste(adj, collapse = '+'), '\n',
                         'z_MUSE_601 ~ b *', pred,  '+ c * lv +', paste(adj, collapse = '+'), '\n',
                         out, '~ d * lv + e * z_MUSE_601 + f *', pred, '+', paste(adj, collapse = '+'))}
  
  if(length(str_dec(ind)) == 1){
    base_model <- paste0(ind, '~ a *', pred, '+', paste(adj, collapse = '+'), '\n',
                         'z_MUSE_601 ~ b *', pred,  '+ c *', ind, '+', paste(adj, collapse = '+'), '\n',
                         out, '~ d *', ind, '+ e * z_MUSE_601 + f *', pred, '+', paste(adj, collapse = '+'))}

    final_model <- paste0(base_model, '\n',
                          '# direct effect
        direct := f
        # indirect effect svd
        in_svd := a * d
        in_svd_gm := a * c * e
        indirect_svd := in_svd + in_svd_gm
        # indirect effect GM
        indirect_gm := b * e
        # total indirect
        indirect_total := indirect_svd + indirect_gm
        # total effect
        total := direct + indirect_total')

  if(c.int == "montCI"){
    fit <- lavaan::sem(final_model, data = data, std.lv = TRUE, estimator = "ML", se = "robust.sem")
    set.seed(2023)
    mat <- semTools::monteCarloCI(object = fit, nRep = rep)}
  
  if(c.int == "boot"){
    fit <- lavaan::sem(final_model, data = data, std.lv = TRUE, estimator = "ML", se = "bootstrap",
                       bootstrap = rep, parallel = "snow", ncpus = (parallel::detectCores() -1), iseed = 2023)
    mat <- lavaan::parameterEstimates(fit, ci = TRUE,boot.ci.type = "bca.simple",
                                      level = 0.95, standardized = TRUE)
    col <- c("lhs", "est", "ci.lower", "ci.upper")
    mat <- mat[mat$op == ":=", colnames(mat) %in% col]
    mat <- mat[,col]
    rownames(mat) <- mat[, "lhs"]
    mat[, c("lhs")] <- NULL}
  
    mat <- mat %>% mutate(across(where(is.numeric), function(x) {
      case_when(abs(x)>=0.001 ~ formatC(x, format = "f", digits = 3),
                TRUE ~ formatC(x, format = "f", digits = 4))}))
    mat$est_ci <- paste0(mat$est, " (", mat$ci.lower, " - ", mat$ci.upper, ")")

  if(r.sq == TRUE){
    rsq <- summary(fit, rsquare = TRUE)$pe
    mat["r2",] <- round(rsq[rsq[,"op"] == "r2" & rsq[,"lhs"] == out, "est"], 3)}
  
  if(sig == TRUE){
    mat$Sig <- ifelse(data.table::between(x = 0,
                                          lower = mat$ci.lower, upper = mat$ci.upper,
                                          incbounds = TRUE), "", "YES")}
  
  mat[, c("est", "ci.lower", "ci.upper")] <- NULL
  colnames(mat)[colnames(mat) == "est_ci"] <- paste(str_dec(ind), collapse = ", ")
  mat <- as.data.frame(t(mat))
  assign(paste(out, paste(str_dec(ind), collapse = ", "), sep = "_"), mat, envir=globalenv())
}

sem_mod_loop <- function(indicators, cor = NULL, out, pred, adj, data, c.int = "montCI",
                         rep = 5000, r.sq = FALSE, sig = FALSE){
  
  str_dec <- function(x) {y <- stringr::str_split_1(x, pattern = "(\\W+\\W+\\W+|\\W+\\W+|\\W+)")
  z <- stringr::str_subset(y, pattern = "(\\w+)")
  return(z)}
  
  nm <- sapply(X = indicators, FUN = function(x) {paste(str_dec(x), collapse = ", ")})
  for (j in seq_along(indicators))
  {cat("> Analysis ", j, " of ", length(indicators), " <", "\nIndicators: ",
       paste(str_dec(indicators[j]), collapse = ", "), "\n\n")
    sem_mod(ind = paste0(indicators[j]), cor = cor, out = out, pred = pred, adj = adj,
            data = data, c.int = c.int, rep = rep, r.sq = r.sq, sig = sig)
    if(all(paste(out, nm, sep = "_") %in% ls(pattern = paste(nm, collapse = "|"), envir = globalenv()))){
      sem_res <- do.call(bind_rows, mget(paste(out, nm, sep = "_"), envir = globalenv()), envir = globalenv())
      assign("sem_res", sem_res, envir = globalenv())
      rm(list = ls(pattern = paste(out), envir = globalenv()), envir = globalenv())
    }}}

ind1 <- "z_bgpvsc"
ind2 <- "z_thpvsc"
ind3 <- "z_pvwml"
ind4 <- "z_otwml"
ind5 <- "z_FA_MUSE_604"
ind6 <- "z_TR_MUSE_604"
ind7 <- "z_bgpvsc + z_thpvsc + z_pvwml + z_otwml + z_FA_MUSE_604 + z_TR_MUSE_604"

# Run function with age as the independent variable
sem_mod_loop(indicators = c(ind1, ind2, ind3, ind4, ind5, ind6, ind7),
             cor = c("z_bgpvsc ~~ z_thpvsc", "z_pvwml ~~ z_otwml", "z_FA_MUSE_604 ~~ z_TR_MUSE_604"),
             out = 'GLOBAL_COG', pred = 'AGE',
             adj = c("SEX_M", "z_MUSE_702", site, race, edu),
             data = mri6, c.int = "boot", rep = 5000, r.sq = TRUE, sig = TRUE)   

# Run function with FRS as the independent variable
sem_mod_loop(indicators = c(ind1, ind2, ind3, ind4, ind5, ind6, ind7),
             cor = c("z_bgpvsc ~~ z_thpvsc", "z_pvwml ~~ z_otwml", "z_FA_MUSE_604 ~~ z_TR_MUSE_604"),
             out = 'GLOBAL_COG', pred = 'frci086c',
             adj = c("z_MUSE_702", site, race, edu),
             data = mri6ct, c.int = "boot", rep = 5000, r.sq = TRUE, sig = TRUE)  

# Multi-group CFA (MGCFA) to test measurement invariance across different sites
# Configural invariance model
model1 <- paste0(
  '# latent variable, measurement model
   svd =~ NA*z_bgpvsc + z_thpvsc + z_pvwml + z_otwml + z_FA_MUSE_604 + z_TR_MUSE_604
   
   # variances
   svd ~~ 1*svd

   # residual covariances
   z_pvwml ~~ z_otwml
   z_bgpvsc ~~ z_thpvsc
   z_FA_MUSE_604 ~~ z_TR_MUSE_604')

# Compute MGCFA model with mean-and-variance corrected chi-square (simple second-order correction) and robust standard errors
fit1 <- lavaan::sem(model1, data = mri6, estimator = "ML", group = "site6c", test = "scaled.shifted", se = "robust.sem")
summary(fit1, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Metric invariance model
model2 <- paste0(
  '# latent variable, measurement model
   svd =~ NA*lbg*z_bgpvsc + lth*z_thpvsc + lpv*z_pvwml + lot*z_otwml + lfa*z_FA_MUSE_604 + ltr*z_TR_MUSE_604
   
   # variances
   svd ~~ c(1, NA, NA, NA, NA, NA)*svd

   # residual covariances
   z_pvwml ~~ z_otwml
   z_bgpvsc ~~ z_thpvsc
   z_FA_MUSE_604 ~~ z_TR_MUSE_604')

# Compute MGCFA model with mean-and-variance corrected chi-square (simple second-order correction) and robust standard errors
fit2 <- lavaan::sem(model2, data = mri6, estimator = "ML", group = "site6c", test = "scaled.shifted", se = "robust.sem")
summary(fit2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Equal inter-site/scanner reliability model
model3 <- paste0(
  '# latent variable, measurement model
   svd =~ NA*lbg*z_bgpvsc + lth*z_thpvsc + lpv*z_pvwml + lot*z_otwml + lfa*z_FA_MUSE_604 + ltr*z_TR_MUSE_604
   
   # variances
   svd ~~ 1*svd
   
   # residual variances
   z_bgpvsc ~~ rbg*z_bgpvsc
   z_thpvsc ~~ rth*z_thpvsc
   z_pvwml ~~ rpv*z_pvwml
   z_otwml ~~ rot*z_otwml
   z_FA_MUSE_604 ~~ rfa*z_FA_MUSE_604
   z_TR_MUSE_604 ~~ rtr*z_TR_MUSE_604

   # residual covariances
   z_pvwml ~~ z_otwml
   z_bgpvsc ~~ z_thpvsc
   z_FA_MUSE_604 ~~ z_TR_MUSE_604')

# Compute MGCFA model with mean-and-variance corrected chi-square (simple second-order correction) and robust standard errors
fit3 <- lavaan::sem(model3, data = mri6, estimator = "ML", group = "site6c", test = "scaled.shifted", se = "robust.sem")
summary(fit3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# Likelihood ratio test with a scaled chi-square difference test statistic
lavaan::lavTestLRT(fit1, fit2, fit3, method = "satorra.2000")

#//----------------------------------------------------------------------------END OF SOURCE CODE----------------------------------------------------------------------------//


