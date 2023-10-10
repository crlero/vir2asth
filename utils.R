if(!require("forestploter", quietly = TRUE))
  install.packages("forestploter")

if(!require("mediation", quietly = TRUE))
  install.packages("mediation")

if(!require("metagMisc", quietly = TRUE))
  remotes::install_github("vmikk/metagMisc")

if(!require("copiome", quietly = TRUE))
  remotes::install_github("jonathanth/copiome@main")

if(!require("tableone", quietly = TRUE))
  devtools::install_github(repo = "kaz-yos/tableone", ref = "develop")

if(!require("mixOmics", quietly = TRUE))
  remotes::install_github("mixOmicsTeam/mixOmics")

if(!require("mixOmicsCaret", quietly = TRUE))
  devtools::install_github("jonathanth/mixOmicsCaret")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("speedyseq", quietly = TRUE))
  BiocManager::install("GenomeInfoDbData", version = '3.15')

if(!require("speedyseq", quietly = TRUE))
  remotes::install_github("mikemc/speedyseq")

if(!require("microViz", quietly = TRUE))
  devtools::install_github("david-barnett/microViz@0.9.4")

if (!require("survMisc", quietly = TRUE))
  install.packages("survMisc")

if (!require("performance", quietly = TRUE))
  install.packages("performance")

packages <- c("tidyverse", # data wrangling
              #"readr", "dplyr", "stringr", "purrr", "tidyr", # data wrangling
              "phyloseq", "vegan", "microbiome", "zCompositions", "metagMisc",
              "copiome", # internal copsac library
              "doParallel", "foreach", "doMC", # parallel computing
              "tableone", # summarise
              "lmerTest", "performance", # statistics
              "ggthemes", "ggsci", "forcats", # visz
              "ggrepel", "ggpubr",
              #"DAtest", # DA benchmarking
              "caret", "mixOmics", "mixOmicsCaret", "pROC") # ML

for (p in packages) {
  library(p, character.only = TRUE)
}

## named vectors
map_urban_to_class <- c("anellovirus" = "Arfiviricetes",
                        "caudovirus" = "Caudoviricetes",
                        "inovirus" = "Faserviricetes",
                        "microvirus" = "Malgrandaviricetes",
                        "adenovirus" = "Tectiliviricetes")

map_class_to_urban <- c("Arfiviricetes" = "anellovirus",
                        "Caudoviricetes" = "caudovirus",
                        "Faserviricetes" = "inovirus",
                        "Malgrandaviricetes" = "microvirus",
                        "Tectiliviricetes" = "adenovirus")

## color vectors
class_color <- c("Microvirus" = "#b16fa9", #"#f7b09a",
                 #"Anellovirus" = "#83cee0",
                 "Caudovirus" = "#52b39f",
                 "Inovirus" = "#f9ba00")
#"Adenovirus" = "#88496f")

endpoint.labs <- c("j45_3yr_cross" = "3 yr (cross-sect.)", # nolint
                   "j45_3yr_ever" = "3 yr (ever)",
                   "j45_3yr_transient" = "3 yr (trans.)",
                   "j45_5yr_cross" = "5 yr (cross-sect.)",
                   "j45_5yr_ever" = "5 yr (ever)",
                   "j45_5yr_transient" = "5 yr (trans.)",
                   "j45_6yr_cross" = "6 yr (cross-sect.)",
                   "j45_6yr_ever" = "6 yr (ever)",
                   "j45_6yr_transient" = "6 yr (trans.)",
                   "oldchild01" = "Siblings",
                   "asthma_mother" = "Asthmatic mother")

dicotom.labs <- c("0" = "Healthy", "1" = "Asthma")
lifestyle.labs <- c(" 0" = "Temperate", " 1" = "Virulent")

taxon.labs <- c("species" = "Species", "genus" = "Genus", "family" = "Family") # nolint

grouping.labs <- c("all" = "Gut virome", "eukaryotic" = "Eukaryotic virome",
                   "bacteriophages" = "Gut bacteriophages",
                   "temperate" = "Temperate virome",
                   "virulent"="Virulent virome")

#color_map <- c("0" = "deepskyblue4", "1" = "brown3")
color_map <- c("0" = "#61a0af", "1" = "#f6511d")
color_rev_map <- c("0" = "#ffffff", "1" = "#ffffff")

color_map2 <- c("Healthy" = "#61a0af",
                "Asthma" = "#f6511d")

color_rev_map2 <- c("Healthy" = "#ffffff",
                    "Asthma" = "#ffffff",
                    "No-siblings" = "#ffffff",
                    "Siblings" = "#ffffff",
                    "No-asthma-mother" = "#ffffff",
                    "Asthma-mother" = "#ffffff")

## functional vectors
lvl <- c("OTU", "category", "class", "phylum", "order", "family",
         "famid", "subfamily", "genus", "species")

## subset phyloseq by lifestyle
subset_lifestyle = function(ps, virclass) {
  if (virclass == "all") {
    print("No viral class subsetting")
    ps
  } else if (virclass == "virulent") {
    print("Subsetting to virulent")
    ps <-  subset_taxa(ps, virulence == " 1")
    ps <- transform_phy(ps, transform = "compositional") # re-normalize
    return(ps)
  } else if (virclass == "temperate") {
    print("Subsetting to temperate")
    ps <-  subset_taxa(ps, virulence == " 0")
    ps <- transform_phy(ps, transform = "compositional") # re-normalize
    return(ps)
  } else if (virclass == "unknown") {
    print("Subsetting to temperate")
    ps <-  subset_taxa(ps, is.na(virulence))
    ps <- transform_phy(ps, transform = "compositional") # re-normalize
    return(ps)
  } else if (virclass == "eukaryote") {
    print("Subsetting to Eukaryome")
    ps <- subset_taxa(ps, class %in% c("Arfiviricetes", "Tectiliviricetes"))
    ps <- transform_phy(ps, transform = "compositional") # re-normalize
  }
}

## calculate distances
calculate_distance = function(ps, threads, midistance, outfile) {
  registerDoParallel(cores = threads)
  if (midistance %in% c("bray", "canberra")) {
    print(paste0("Distance: ", midistance))
    dist <- distance(ps, method = midistance, parallel = TRUE)
    saveRDS(list(dist, ps), file = outfile)
  } else if (midistance %in% c("jaccard", "unifrac")) {
    print(paste0("Distance: ", midistance))
    ps <- transform_phy(ps, transform = "pa")
    if (midistance == "jaccard") {
      dist <- distance(ps, method = midistance,
                       binary = TRUE, parallel = TRUE)
    } else {
      dist <- distance(ps, method = midistance, parallel = TRUE)
    }
    saveRDS(list(dist, ps), file = outfile)
  } else if (midistance == "euclidean") {
    print(paste0("Distance: ", midistance))
    ps_clr <- transform_phy(ps, transform = "clr")
    dist <- distance(ps_clr, method = midistance, parallel = TRUE)
    saveRDS(list(dist, ps_clr), file = outfile)
  } else if (midistance %in% c("jsd", "wunifrac")) {
    print(paste0("Distance: ", midistance))
    ps <- transform_phy(ps, transform = "log")
    dist <- distance(ps, method = midistance, parallel = TRUE)
    if (midistance == "jsd") {
      dist <- sqrt(dist)
    }
    saveRDS(list(dist, ps), file = outfile)
  } else if (midistance == "combined") {
    print(paste0("Distance: ", midistance))
    dist_bray <- distance(ps, method = "bray",
                          parallel = TRUE) # bray-curtis
    dist_unifrac <-  distance(transform_phy(ps, transform = "pa"),
                              method = "unifrac",
                              parallel = TRUE)
    dist <- combined_metric(dist_bray, dist_unifrac, alpha = 0.5)
    saveRDS(list(dist, ps), file = outfile)
  } else if (midistance == "unifrac-jsd") {
    print(paste0("Distance: ", midistance))
    dist_unifrac <- distance(ps, method = "unifrac", parallel = TRUE)
    ps_log <- transform_phy(ps, transform = "log")
    dist_jsd <- sqrt(distance(ps_log, method = "jsd", parallel = TRUE))
    dist <- dist_unifrac * dist_jsd
    saveRDS(list(dist, ps_log), file = outfile)
  }
}


## calculate alpha diversity metrics
calculate_alpha <- function(ps) {
  if (nrow(tax_table(ps)) > 1) {
    print("Calculating alpha diversity indices...")
    otu_table(ps) <- otu_table(round(otu_table(ps)*1e9)) # use absolute counts, round them (i.e. for shannon and chao)
    df <- estimate_richness(ps,
                            measures = c("Observed", "Shannon", "Chao1")) %>%
      mutate(abcno = sample_df(ps)$abcno) %>% 
      left_join(sample_df(ps), .)
    return(df)
  } else {
    print("Only one tip in the tree.
          No alpha diversity will be calculated at this level.")
  }
}

extract_alpha <- function(ps, group) {
  if (group == "all") {
    print("No subsetting")
    print("Calculating alpha diversity")
    df_richness <- calculate_alpha(ps) %>%
      rename_with(.cols = c("Observed", "Shannon", "Chao1", "se.chao1"),
                  ~paste0(group, "_", tolower(.x)))
    return(df_richness)
  } else if (group == "virulent") {
    print("Subsetting to virulent")
    ps_tmp <-  subset_taxa(ps, virulence == " 1") %>%
      transform(., transform = "compositional")
    print("Calculating alpha diversity")
    df_richness <- calculate_alpha(ps_tmp) %>%
      rename_with(.cols = c("Observed", "Shannon", "Chao1", "se.chao1"),
                  ~paste0(group, "_", tolower(.x)))
    return(df_richness)
  } else if (group == "temperate") {
    print("Subsetting to temperate")
    ps_tmp <-  subset_taxa(ps, virulence == " 0") %>%
      transform(., transform = "compositional")
    print("calculating alpha diversity")
    df_richness <- calculate_alpha(ps_tmp) %>%
      rename_with(.cols = c("Observed", "Shannon", "Chao1", "se.chao1"),
                  ~paste0(group, "_", tolower(.x)))
    return(df_richness)
  } else if (group == "unknown") {
    print("Subsetting to uknown lifestyle")
    ps_tmp <- subset_taxa(ps, is.na(virulence)) %>% 
      transform(., transform = "compositional")
    print("calculating alpha diversity")
    df_richness <- calculate_alpha(ps_tmp) %>% 
      rename_with(.cols = c("Observed", "Shannon", "Chao1", "se.chao1"),
                  ~paste0(group, "_", tolower(.x)))
  }
}

## associations alpha - outcome
associate_alpha_outcome <- function(tmp, vars, div, seed=123) {
  mp_list <- list()
  emm_list <- list()
  for (var in vars) {
    print(paste0("Analyzing: ", var, " - ", div))
    tmp_t <- tmp %>%
      mutate_at(.vars = c(var), as.factor) %>%
      filter_at(var, all_vars(!is.na(.)))
    if (grepl("transient", var)) {
      tmp_t <- tmp_t[tmp_t[[stringr::str_replace(var, "transient", "ever")]] == 1 | tmp_t[[var]] == 1, ] # no-lint
    }
    set.seed(seed) # regression
    f <- as.formula(paste(div, "~",
                          "+ log(mapDepth) + mapDepth + efficiency \ 
                          + propOTU + age2010 + lane"))
    fit <- lm(f, tmp_t)
    tmp_t$calibrated <- coef(fit)[1] + resid(fit)
    
    res_unadj <- wilcox.test(as.formula(paste("calibrated ~ ", var)),
                             data = tmp_t)
    mp <- model_parameters(res_unadj) %>%
      data.frame() %>%
      mutate(diversity = div,
             outcome = var, 
             N = nrow(tmp_t),
             n = nrow(tmp_t[tmp_t[[var]] == 1, ]))
    mp_list[[var]] <- mp
    # calibrate
    M <- tmp_t %>% group_by_at(.vars = var) %>%
      summarise_at(.funs = mean, .vars = "calibrated")
    M$mean <- paste0("M=", round(M[["calibrated"]], 2))
    # median
    Mdn <- tmp_t %>% group_by_at(.vars = var) %>% summarise_at(.funs = median, .vars = "calibrated")
    Mdn$median <- paste0("Mdn=", round(Mdn[["calibrated"]], 2))
    # stats
    summ <- M %>% dplyr::select_at(., .vars = c(var, "mean")) %>%
      left_join(., dplyr::select_at(Mdn, .vars = c(var, "median"))) %>%
      mutate(label = paste0(mean, "\n", median))
    # figure
    figboxplot <- ggplot(tmp_t, aes_string(x = var, y = "calibrated", fill = var)) +
      geom_point(position = position_jitter(w = 0.1, h = 0), pch = 21,
                 size = 2.5, aes_string(fill = var, color = var)) +
      geom_boxplot(alpha = 0.5, outlier.size = -Inf,
                   width = 0.5 / length(unique(tmp_t[[var]]))) +
      theme_minimal(base_size = 12) +
      xlab("") + ylab(div) +
      guides(fill = "none", color = "none") +
      scale_color_manual(values = color_rev_map) +
      scale_fill_manual(values = color_map) +
      scale_x_discrete(labels = c("0" = "Healthy", "1" = "Asthma")) +
      geom_label(aes_string(label = "label"), y = 1.05 * max(tmp_t[["calibrated"]]),
                 color = "white", data = summ, size = 3) +
      expand_limits(y = 1.05 * max(tmp_t[["calibrated"]])) +
      labs(title = paste0(endpoint.labs[var], " (", nrow(tmp_t[tmp_t[[var]] == 1,]),
                          "/", nrow(tmp_t), ")"),
           subtitle = paste0("Wilcoxon, P = ",
                             formatC(res_unadj$p.value, format = "e", digits = 1)))
    emm_list[[var]] <- figboxplot
    if (res_unadj$p.value < 0.05) {
      print(paste0(var, " ~ ", div, " is significant"))
    }
  }
  mp_df <- bind_rows(mp_list) %>%
    arrange(desc(-p))
  batchvars <- c("lane", "log(mapQcDepth)",
                 "efficiency", "log(mapDepth)",
                 "propViral", "log(qcDepth)", "propOTU", "viromeQC")
  return(list(mp_df, emm_list))
}

## Remove samples when phenotype is NA
remove_na_samples <- function(ps, var_name) {
  var_values <- phyloseq::sample_data(ps)[[var_name]]
  phyloseq::prune_samples(!is.na(var_values), ps)
}


## Plot prevalence vs. mra
prevalence_mra_plot <- function(ps) {
  mra <- make_mradat(ps) %>%
    left_join(., prevalence_df(ps)) %>% 
    mutate(class_urban = str_to_sentence(map_class_to_urban[as.character(class)]))
  cor.res <- cor.test(mra$prevalence, mra$mra, method="spearman", exact=FALSE)
  cor.res
  miorder <- mra %>% group_by(class_urban) %>% summarise(top=n()) %>% arrange(desc(top)) %>% .$class_urban
  miorder <- data.frame(class_urban=miorder, miorder=seq(1:length(miorder)))
  miorder <- miorder %>% mutate(miorder = ifelse(class_urban == "Microvirus", 2, miorder),
                                miorder = ifelse(class_urban == "Inovirus", 3, miorder))
  mra <- mra %>% left_join(., miorder) %>% arrange(miorder, -mra)
  mra <- mra %>% mutate(class = as.character(class),
                        class_urban = str_to_sentence(map_class_to_urban[class]),
                        class_urban = factor(class_urban, levels=c("Caudovirus", "Microvirus", "Inovirus")))
  
  p <- ggplot(mra, aes(x=mra*100, y=prevalence)) +
    geom_point(size=2, alpha=.5, pch=21, #color="white", lwd=0.05,
               stroke=.2, aes(fill=class_urban, color=class_urban), data = . %>% filter(class_urban == "Caudovirus")) +
    geom_point(size=2, alpha=.5, pch=21, #color="white", lwd=.05,
               stroke=.2, aes(fill=class_urban, color=class_urban), data = . %>% filter(class_urban == "Microvirus")) +
    geom_point(size=2, alpha=.5, pch=21, #color="white", lwd=.05,
               stroke=.2, aes(fill=class_urban, color=class_urban), data = . %>% filter(class_urban == "Inovirus")) +
    #geom_rug(aes(color=class_urban)) +
    scale_x_continuous(trans='log10', labels = scales::comma) +
    scale_fill_manual(values=class_color, name="Class") +
    scale_color_manual(values=class_color, name="Class") +
    scale_size_continuous(name="N vOTUs") +
    # guides(fill="none") +
    theme_minimal() +
    xlab("mra (%)") +
    ylab("Prevalence (%)")
  return(p)
}


mydbrda <- function(dist, ps, seed=123, verbose=TRUE) {
  phenodata <- sample_df(ps)
  set.seed(seed)
  rdares <- dbrda(dist ~ 1 + Condition(lane), phenodata)
  if (verbose) {
    print(rdares)
  }
  return(rdares)
}


plot_beta <- function(dist, ps, var, nproc=4, distance_metric="", g="",
                      verbose=TRUE, formula=NA, seed=123) {
  phenodata <- sample_df(ps) %>%
    filter_at(., vars(var), all_vars(!is.na(.)))
  res_strata_df <- myadonis(dist, ps, var, nproc=nproc)
  rdares <- mydbrda(dist, ps, verbose = TRUE, seed = seed)
  axis1 <- round(rdares$CA$eig[["MDS1"]] / rdares$tot.chi * 100, 2)
  axis2 <- round(rdares$CA$eig[["MDS2"]] / rdares$tot.chi * 100, 2)
  scores <- data.frame(scores(rdares, display = "sites")) %>%
    tibble::rownames_to_column(., "abcno")
  df <- phenodata %>%
    left_join(., scores)
  mdsplot.all <- ggplot(df[!is.na(df[[var]]), ],
                        aes_string(x = "MDS1", y = "MDS2",
                                   fill = var, color = var)) +
    stat_ellipse(aes_string(color = var, fill = var), geom = "polygon", alpha=.2, size=.20, type="t") +
    geom_point(pch = 21, stroke = .3, size = 2.5, alpha = .9) +
    scale_fill_manual(values = color_map) + scale_color_manual(values =  color_rev_map) +
    xlab(paste0("MDS1 (", axis1, "%)")) +
    ylab(paste0("MDS2 (", axis2, "%)")) +
    #labs(title = var,
    labs(title = g,
         caption = distance_metric,
         subtitle = paste0("PERMANOVA, R2(%): ",
                           round(subset(res_strata_df,
                                        parameter == var)$R2 * 100, 2),
                           "; P = ",
                           format(subset(res_strata_df,
                                         parameter == var)$pval,
                                  format = "e", digits = 1))) +
    guides(fill="none", color="none") +
    theme_bw()
  return(mdsplot.all)
}


get_loadings = function (trainobj, what = c("finalModel", "CV"), xykeep = c("x", 
                                                                            "y", "both"), rep = NA, remove_empty = TRUE, summarize = TRUE, 
                         ncomp = NA, keepX = NA, keepY = NA) 
{
  if (what[1] == "finalModel") {
    loadingsx <- trainobj$finalModel$loadings$X %>% data.frame %>% 
      mutate(var = rownames(.), sd = NA, chosen = NA) %>% 
      gather(comp, loading, -var, -sd)
    loadingsy <- trainobj$finalModel$loadings$Y %>% data.frame %>% 
      mutate(var = rownames(.), sd = NA, chosen = NA) %>% 
      gather(comp, loading, -var, -sd)
  }
  if (what[1] == "CV") {
    bt <- trainobj$bestTune
    if (!is.na(ncomp)) 
      bt$ncomp <- ncomp
    if (!is.na(keepX)) 
      bt$keepX <- keepX
    if (!is.na(keepY)) 
      bt$keepY <- keepY
    control <- trainobj$control
    keep <- rep(TRUE, length(control$index))
    if (!is.na(rep)) 
      keep <- grepl(rep, names(control$index))
    index <- control$index[keep]
    indexOut <- control$indexOut[keep]
    seeds <- control$seeds[keep]
    refits <- lapply(seq(along = index), function(i) {
      trainobj$modelInfo$fit(trainobj$finalModel$X[index[[i]],],
                             trainobj$finalModel$Y[index[[i]], ],
                             param = bt, 
                             trainobj$dots)
    })
    loadingsx <- lapply(seq_along(refits), function(i) {
      refits[[i]]$loadings$X %>% data.frame %>% mutate(var = rownames(.)) %>% 
        gather(comp, loading, -var)
    }) %>% bind_rows
    loadingsy <- lapply(seq_along(refits), function(i) {
      refits[[i]]$loadings$Y %>% data.frame %>% mutate(var = rownames(.)) %>% 
        gather(comp, loading, -var)
    }) %>% bind_rows
    if (summarize) {
      loadingsx <- loadingsx %>% group_by(comp, var) %>% 
        dplyr::summarize(sd = sd(loading), loading = median(loading), 
                         chosen = mean(loading != 0)) %>% ungroup
      loadingsy <- loadingsy %>% group_by(comp, var) %>% 
        dplyr::summarize(sd = sd(loading), loading = median(loading), 
                         chosen = mean(loading != 0)) %>% ungroup
    }
  }
  if (!what[1] %in% c("finalModel", "CV")) {
    stop("Please choose \"finalModel\" or \"CV\".")
  }
  out <- dplyr::bind_rows(x = loadingsx, y = loadingsy, .id = "xy") %>% 
    dplyr::select(var, comp, loading, sd, xy)
  if (remove_empty) 
    out <- out %>% dplyr::filter(loading != 0)
  if (xykeep[1] %in% c("x", "y")) 
    return(out %>% filter(xy == xykeep[1]))
  out
}

## Plot PLS-DA model components AUC distribution
auc_components_plot = function(plsmodel) {
  plsmodel$pred %>%
    separate(Resample, c("Fold", "Rep")) %>%
    group_by(ncomp, keepX, Rep) %>%
    summarize(auc = as.numeric(pROC::auc(predictor = pred, obs, direction = "<"))) %T>%
    { mx <<- max(.$auc); mn <<- min(.$auc) } %>%
    ggplot(., aes(x = factor(keepX), y = auc, fill=factor(ncomp))) +
    geom_boxplot(outlier.shape = NA, alpha=0.25) +
    geom_point(aes(color=factor(ncomp)),
               alpha=0.6,
               position=position_jitter(w=0.15, h=0)) +
    guides(fill="none", color="none") +
    facet_wrap(~ ncomp) +
    geom_hline(yintercept = 0.5, lwd=.5, linetype=2) +
    scale_color_brewer(palette = "Set1", name = NULL) +
    theme_bw() + theme(strip.background = element_blank()) +
    xlab("Number of features") + ylab("AUC")
}

plot_procrustes = function(procrustes_result, nround = 2, var = NULL) {
  pro.dat <- data.frame(mds1 = procrustes_result$Yrot[, 1],
                        mds2 = procrustes_result$Yrot[, 2],
                        xmds1 = procrustes_result$X[, 1],
                        xmds2 = procrustes_result$X[, 2])
  if (!is.null(var)) {
    pro.dat$group <- var
  }
  pl <- ggplot(pro.dat)
  if (!is.null(var)) {
    pl <-  pl + 
      geom_segment(aes(x=xmds1, y=xmds2, xend=mds1, yend=mds2, color = as.factor(group)),
                   size=0.25,
                   arrow=arrow(length=unit(0.15,"cm"), ends="first"), show.legend=F) +
      geom_point(aes(x=mds1, y=mds2, color = as.factor(group)), size=1, pch=21) +
      scale_color_manual(values = color_map)
  } else {
    pl <-  pl + 
      geom_segment(aes(x=xmds1, y=xmds2, xend=mds1, yend=mds2),
                   size=0.25, color="steelblue",
                   arrow=arrow(length=unit(0.15,"cm"), ends="first"), show.legend=F) +
      geom_point(aes(x=mds1, y=mds2), color = "steelblue", size=1, pch=21)
  }
  pl <- pl + 
    grids(linetype = "dashed") +
    geom_hline(yintercept = 0, linetype="dashed", lwd = .25) +
    geom_vline(xintercept = 0, linetype="dashed", lwd = .25) +
    geom_abline(slope = -1/procrustes_result$rotation[1,2], lwd = .25) +
    geom_abline(slope = procrustes_result$rotation[1,2], lwd = .25) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.05)) + # to adjust decimals
    scale_y_continuous(labels = scales::number_format(accuracy = 0.05)) +
    theme(panel.background = element_blank(),
          legend.position = "bottom",
          legend.background = element_blank(),
          legend.key = element_blank(),
          axis.line = element_line()) +
    xlab("MDS1") + ylab("MDS2") +
    annotate("text", Inf, -Inf, label = paste0("P == ", procrustes_result$signif),
             parse = TRUE, size = 4, hjust = 1.1, vjust = -1.6) +
    annotate("text", Inf, -Inf, label = paste0("italic(r) == ", roundex(procrustes_result$scale, nround)),
             parse = TRUE, size = 4, hjust = 1.1, vjust = -0.4) #+
  # annotate("text", Inf, -Inf, label = paste0("italic(m) ^ 2 == ", roundex(procrustes_result$ss, nround)),
  #          parse = TRUE, size = 4, hjust = 1.1, vjust = -3.8)
  return(pl)
}


# reorder samples
ps_reorder <- function(ps, sample_order) {
  otu <- phyloseq::otu_table(ps)
  if (phyloseq::taxa_are_rows(ps)) {
    otu <- otu[, sample_order]
  } else {
    otu <- otu[sample_order, ]
  }
  phyloseq::otu_table(ps) <- otu
  
  return(ps)
}


signif_stars <- function(x) {
  unclass(
    symnum(x,
           corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
           symbols = c("***", "**", "*", ".", " ")
    )
  )
}

# extract files
get_midf <- function(dir, ncores=3, resolution, group) {
  registerDoParallel(cores=ncores)
  df <- foreach(f=dir, .combine=bind_rows) %dopar% {
    dist_name <- str_replace(tail(unlist(str_split(f, "/")), 1), ".tsv", "")
    read_tsv(f) |> mutate(resolution = resolution, 
                          group = group, distance = dist_name,
                          significant = signif_stars(pval)) |> 
      filter(variable == "j45_5yr_ever")
  }
  return(df |> arrange(pval))
}


# nice visz for tables
dt <- function(d, ...) {
  datatable(d, 
            extensions = "FixedColumns",
            filter = list(position = "top", clear = TRUE), 
            rownames = FALSE, 
            options = c(list(pageLength = 5, autoWidth = TRUE, scrollX = TRUE, 
                             searching = TRUE, lengthChange = FALSE, 
                             search = list(regex = TRUE),
                             columnDefs = list(list(width = "120px", targets = "_all")),
                             fixedColumns = list(leftColumns = 1)), ...))
}
