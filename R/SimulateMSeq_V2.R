##############################################
########## Simulate Data Functions ###########
##############################################
## Dirichlet Function
rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}


## Estimate the Dirichlet Hyper-parameters Using MLE
EstPara <- function (ref.otu.tab) {
  
  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0('OTU', 1 : nrow(ref.otu.tab))
  } # otu * sample
  samplenames = colnames(ref.otu.tab)
  taxnames = rownames(ref.otu.tab)
  
  dirmult.paras <- dirmult::dirmult(t(ref.otu.tab))
  
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  ref.otu.tab = sapply(1:ncol(ref.otu.tab), function (i) gamma + ref.otu.tab[,i]) # C_ij otu * sample
  
  # back to Dirichlet, calculate the true proportion
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
  colnames(ref.otu.tab.p) = samplenames
  rownames(ref.otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p =  ref.otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 = t(t(ref.otu.tab.p)*Si)
  colnames(ref.otu.tab0) = colnames(ref.otu.tab.p)
  return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}

## Simulate Data with estimated Hyper-parameters
SimulateMSeqU<-function (para = EstPara, nSam = 100, nOTU = 500, diff.otu.pct = 0.1, 
                         diff.otu.direct = c("balanced", "unbalanced"), 
                         diff.otu.mode = c("abundant",  "rare", "mix","user_specified"), 
                         user_specified_otu = NULL,
                         covariate.type = c("binary", "continuous"), 
                         grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0, 
                         confounder.type = c("none", "binary", "continuous", "both"), 
                         conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1, 
                         confounder.eff.mean = 0, confounder.eff.sd = 0, error.sd = 0, 
                         depth.mu = 10000, depth.theta = 5, depth.conf.factor = 0, cont.conf, epsilon)
{
  ## Estimated Dirichlet Parameters
  model.paras <- para
  
  ## Select number of OTUs
  ref.otu.tab <- model.paras$ref.otu.tab[(1:(nOTU)), ]
  
  ## --- Select / expand number of Samples ---
  expand_samples <- function(tab, nSam) {
    tab <- as.matrix(tab)
    n_ref <- ncol(tab)
    
    idx <- sample(seq_len(n_ref), size = nSam, replace = TRUE)
    
    out <- tab[, idx, drop = FALSE]
    colnames(out) <- paste0(colnames(tab)[idx], ".sim", seq_len(nSam))
    
    out
  }
  
  ref.otu.tab <- expand_samples(ref.otu.tab, nSam)
  idx.otu     <- rownames(ref.otu.tab)
  idx.sample  <- colnames(ref.otu.tab)
  ## Select number of Samples
  # idx.otu <- rownames(ref.otu.tab)
  # idx.sample <- colnames(model.paras$ref.otu.tab)[1:nSam]
  ref.otu.tab = ref.otu.tab[, idx.sample]
  
  ## Confounding variable
  if (confounder.type == "none") {
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
    Z <- NULL
  }
  if (confounder.type == "continuous") {
    Z <- cbind(cont.conf)
  }
  
  if (confounder.type == "binary") {
    Z <- cbind(c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  }
  
  if (confounder.type == "both") {
    Z <- cbind(rnorm(nSam), c(rep(0, nSam%/%2), rep(1, nSam - nSam%/%2)))
  } 
  
  ## Covariate of Interest
  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  
  if (covariate.type == "continuous") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + epsilon#rnorm(nSam)
  }
  
  if (covariate.type == "binary") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + epsilon#rnorm(nSam)
    X <- cbind(ifelse(X <= quantile(X, grp.ratio/(1 + grp.ratio)), 0, 1))
  }
  
  rownames(X) <- colnames(ref.otu.tab)
  covariate.eff.mean1 = covariate.eff.mean
  covariate.eff.mean2 = covariate.eff.mean
  
  ## Simulate OTU Absolute Abundance
  ## Balanced OTU Counts
  if (diff.otu.direct == "balanced") {
    if (diff.otu.mode == "user_specified") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    if (diff.otu.mode == "abundant") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*% 
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/4), mean = -covariate.eff.mean1, sd = covariate.eff.sd), 
                             rnorm(floor(nOTU/2) - floor(nOTU/4), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(floor((nOTU - floor(nOTU/2))/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                             rnorm(nOTU - floor(nOTU/2) - floor((nOTU - floor(nOTU/2))/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% 
        t(scale(X))
    }
  }
  
  ## Unbalanced OTU Counts
  if (diff.otu.direct == "unbalanced") {
    if (diff.otu.mode == "user_specified") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% 
        t(scale(X))
    }
    if (diff.otu.mode == "abundant") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% 
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% 
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/2), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*%
        t(scale(X))
    }
  }
  eta.conf <- sample(c(rnorm(floor(nOTU/2), mean = -confounder.eff.mean, sd = confounder.eff.sd), 
                       rnorm(nOTU - floor(nOTU/2), mean = confounder.eff.mean, sd = confounder.eff.sd))) %*% 
    t(scale(scale(Z) %*% rep(1, ncol(Z))))
  
  otu.ord <- 1:(nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  
  if (diff.otu.mode == "user_specified") 
    diff.otu.ind <- c(diff.otu.ind, which(idx.otu %in% user_specified_otu))
  if (diff.otu.mode == "mix") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
  if (diff.otu.mode == "abundant") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1:round(length(otu.ord)/4)], diff.otu.num))
  if (diff.otu.mode == "rare") 
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * length(otu.ord)/4):length(otu.ord)], diff.otu.num))
  
  if (length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
    conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
  } else {
    conf.otu.ind1 <- diff.otu.ind
  }
  conf.otu.ind <- c(conf.otu.ind1, 
                    sample(setdiff(1:(nOTU), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
  
  ## Calculate the new composition
  eta.diff[setdiff(1:(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1:(nOTU), conf.otu.ind), ] <- 0
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, nSam)
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  ref.otu.tab.prop <- eta.exp/rowSums(eta.exp) ## Pij'=  Normalized Cij''
  #ref.otu.tab.prop = eta.exp
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  #colnames(ref.otu.tab.prop) <- rownames(eta.exp)
  #rownames(ref.otu.tab.prop) <- rownames(ref.otu.tab)
  
  ## Simulate Sequencing Depth using Negative Binomial Distribution
  nSeq <- rnegbin(nSam, 
                  mu = depth.mu * exp(scale(X) * depth.conf.factor), 
                  theta = depth.theta)
  
  ## Simulate Absolute Abundance of OTU
  otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop), 
                        function(i) rmultinom(1, nSeq[i], ref.otu.tab.prop[, i]))
  #otu.tab.sim <- round(t(t(ref.otu.tab.prop)*nSeq))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  diff.otu.ind = (1:nOTU) %in% diff.otu.ind
  conf.otu.ind = (1:nOTU) %in% conf.otu.ind
  
  ## Output
  return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z, 
              diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}

##############################################
########## Check results Functions ###########
##############################################

## Global Test Function
Global_Test = function(sim.dat){
  OTU = t(sim.dat$otu.tab.sim)
  Group = sim.dat$covariate
  data = data.frame(OTU,Group)
  
  ## Shannon Index and Wilcoxon Test
  data$shannon = vegan::diversity(OTU, index = "shannon")
  shannon_group = split(data$shannon, data$Group)
  Shannon_p = wilcox.test(shannon_group$`0`, shannon_group$`1`, alternative = "two.sided")$p.value
  
  ## Inverse Simpson Index and Wilcoxon Test
  data$invsimpson = vegan::diversity(OTU, index = "invsimpson")
  invsimpson_group = split(data$invsimpson, data$Group)
  Inv_Simpson_P = wilcox.test(invsimpson_group$`0`, invsimpson_group$`1`, alternative = "two.sided")$p.value
  
  ## PERMANOVA Test on Bray Curtis Distances
  bray_curtis_distances = as.matrix(vegdist(OTU, method = "bray"))
  permanova.bray = adonis2(bray_curtis_distances ~ Group, permutations = 999,data = data)
  Bray_Curtis_P = permanova.bray$`Pr(>F)`[1]
  
  out = c(Shannon_p,Inv_Simpson_P,Bray_Curtis_P)
  return(out)
}

## MiRKAT Integration
MiRKAT_Test = function(sim.index){
  
  ############# 16s Kernel Matrix #############
  sim.dat.16s = sim_16S_Genus[[sim.index]]
  OTU.16s = t(sim.dat.16s$otu.tab.sim)
  Group.16s = as.numeric(sim.dat.16s$covariate)
  ## Bray Curtis Distance
  bc.16s = as.matrix(vegdist(OTU.16s, method = "bray"))
  ## converts distance matrices to kernel matrices
  K.BC.16s = D2K(bc.16s)
  
  ############# WGS Kernel Matrix #############
  sim.dat.wgs = sim_WGS_Genus[[sim.index]]
  OTU.wgs = t(sim.dat.wgs$otu.tab.sim)
  Group.wgs = as.numeric(sim.dat.wgs$covariate)
  ## Bray Curtis Distance
  bc.wgs = as.matrix(vegdist(OTU.wgs, method = "bray"))
  ## converts distance matrices to kernel matrices
  K.BC.wgs = D2K(bc.wgs)
  
  
  ## Testing Using a Single Kernel
  Ks = list(K.BC.16s,K.BC.wgs)
  
  MiRKAT_result = MiRKAT(y = Group.wgs, X = NULL, Ks = Ks, out_type = "D",
                         omnibus = "permutation", nperm = 999, 
                         returnKRV = TRUE, returnR2 = TRUE)
  
  
  out = c(MiRKAT_16s_p = ifelse(MiRKAT_result$p_values[1]<0.05,1,0),
          MiRKAT_wgs_p = ifelse(MiRKAT_result$p_values[2]<0.05,1,0),
          MiRKAT_omnibus_p = ifelse(MiRKAT_result$omnibus_p<0.05,1,0))
}

## Univariate Test Function
TP_Rate = function(simulated_dat){
  
  count = simulated_dat$otu.tab.sim
  col = data.frame(Group = factor(simulated_dat$covariate),row.names = rownames(simulated_dat$covariate))
  
  ## Total True Positive features
  Total_P = simulated_dat$otu.names[simulated_dat$diff.otu.ind]
  
  ############################
  ###### DEseq2 Analysis #####
  ############################
  dds = DESeqDataSetFromMatrix(countData = count, colData = col, design = ~ Group)
  DeSeq_results = DESeq(dds, test="Wald",fitType = "parametric",sfType = "poscounts")
  DeSeq_results_df = results(DeSeq_results,pAdjustMethod = "BH")
  
  ## Select Significant 
  select_sig =  !is.na(DeSeq_results_df$padj) & DeSeq_results_df$padj < 0.05
  select_sig_df = DeSeq_results_df[select_sig,]
  select_sig_names = rownames(select_sig_df)
  
  ## Calculate TPR, FDR, Matthews correlation coefficient
  classification = lapply(rownames(count), function(feature){
    actual = feature %in% Total_P
    pred = feature %in% select_sig_names
    return(c(feature, actual, pred))
  })
  
  classification_df = data.frame(do.call(rbind,classification))
  colnames(classification_df) = c("taxon","actual","pred")
  
  DEseq2_TPR = sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(Total_P)
  DEseq2_TNR = sum(classification_df$actual == FALSE & classification_df$pred==FALSE)/(nrow(classification_df) - length(Total_P))
  DEseq2_FDR = 1- sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(select_sig_names)
  DEseq2_MCC = mcc(preds = classification_df$pred,actuals = classification_df$actual)
  
  ############################
  ###### ANCOMBC Analysis ####
  ############################
  # phylo = phyloseq(otu_table(count,taxa_are_rows = TRUE),sample_data(col))
  # anc_results = ancombc(data = phylo,formula = "Group",p_adj_method = "BH")
  # anc_results_df = anc_results[["res"]][["p_val"]]
  # 
  # ## Select Significant 
  # select_sig =  anc_results_df$Group1 < 0.05
  # select_sig_df = anc_results_df[select_sig,]
  # select_sig_names = select_sig_df$taxon
  # 
  # ## Calculate TPR, FDR, Matthews correlation coefficient
  # classification = lapply(rownames(count), function(feature){
  #   actual = feature %in% Total_P
  #   pred = feature %in% select_sig_names
  #   return(c(feature, actual, pred))
  # })
  # 
  # classification_df = data.frame(do.call(rbind,classification))
  # colnames(classification_df) = c("taxon","actual","pred")
  # 
  # ANCOMBC_TPR = sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(Total_P)
  # ANCOMBC_TNR = sum(classification_df$actual == FALSE & classification_df$pred==FALSE)/(nrow(classification_df) - length(Total_P))
  # ANCOMBC_FDR = 1- sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(select_sig_names)
  # ANCOMBC_MCC = mcc(preds = classification_df$pred,actuals = classification_df$actual)
  
  
  ############################
  ###### ANCOMBC2 Analysis ###
  ############################
  # phylo = phyloseq(otu_table(count,taxa_are_rows = TRUE),sample_data(col))
  # anc2_results = ancombc2(data = phylo,fix_formula = "Group",p_adj_method = "BH",pseudo_sens = FALSE)
  # anc2_results_df = anc2_results[["res"]]
  # 
  # ## Select Significant 
  # select_sig =  anc2_results_df$p_Group1 < 0.05
  # select_sig_df = anc2_results_df[select_sig,]
  # select_sig_names = select_sig_df$taxon
  # 
  # ## Calculate TPR, FDR, Matthews correlation coefficient
  # classification = lapply(rownames(count), function(feature){
  #   actual = feature %in% Total_P
  #   pred = feature %in% select_sig_names
  #   return(c(feature, actual, pred))
  # })
  # 
  # classification_df = data.frame(do.call(rbind,classification))
  # colnames(classification_df) = c("taxon","actual","pred")
  # 
  # ANCOMBC2_TPR = sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(Total_P)
  # ANCOMBC2_TNR = sum(classification_df$actual == FALSE & classification_df$pred==FALSE)/(nrow(classification_df) - length(Total_P))
  # ANCOMBC2_FDR = 1- sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(select_sig_names)
  # ANCOMBC2_MCC = mcc(preds = classification_df$pred,actuals = classification_df$actual)
  
  
  ############################
  ##### Wilcoxon Analysis ####
  ############################
  # count1 = data.frame(count[,rownames(col)[col$Group==1]])
  # count2 = data.frame(count[,rownames(col)[col$Group==0]])
  # 
  # wilcoxon_results = lapply(rownames(count),function(feature){
  #   wil = wilcox.test(as.numeric(count1[feature,]),as.numeric(count2[feature,]))
  #   return(c(feature,wil$statistic,wil$p.value))
  # })
  # 
  # wilcoxon_results_df = data.frame(do.call(rbind, wilcoxon_results))
  # colnames(wilcoxon_results_df) = c("taxon","W","p_val")
  # wilcoxon_results_df$p_adj = p.adjust(wilcoxon_results_df$p_val,"BH")
  # 
  # ## Select Significant 
  # select_sig =   !is.na(wilcoxon_results_df$p_adj) & wilcoxon_results_df$p_adj < 0.05
  # select_sig_df =  wilcoxon_results_df[select_sig,]
  # select_sig_names = select_sig_df$taxon
  # 
  # ## Calculate TPR, FDR, Matthews correlation coefficient
  # classification = lapply(rownames(count), function(feature){
  #   actual = feature %in% Total_P
  #   pred = feature %in% select_sig_names
  #   return(c(feature, actual, pred))
  # })
  # 
  # classification_df = data.frame(do.call(rbind,classification))
  # colnames(classification_df) = c("taxon","actual","pred")
  # 
  # wilcoxon_TPR = sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(Total_P)
  # wilcoxon_TNR = sum(classification_df$actual == FALSE & classification_df$pred==FALSE)/(nrow(classification_df) - length(Total_P))
  # wilcoxon_FDR = 1- sum(classification_df$actual == TRUE & classification_df$pred==TRUE)/length(select_sig_names)
  # wilcoxon_MCC = mcc(preds = classification_df$pred,actuals = classification_df$actual)
  
  ############################
  ####### Return Results #####
  ############################
  
  # return(c(DEseq2_TPR, DEseq2_TNR, DEseq2_FDR, DEseq2_MCC,
  #          ANCOMBC_TPR, ANCOMBC_TNR, ANCOMBC_FDR, ANCOMBC_MCC,
  #          ANCOMBC2_TPR, ANCOMBC2_TNR, ANCOMBC2_FDR,ANCOMBC2_MCC,
  #          wilcoxon_TPR, wilcoxon_TNR, wilcoxon_FDR,wilcoxon_MCC))
  return(c(DEseq2_TPR))
}


