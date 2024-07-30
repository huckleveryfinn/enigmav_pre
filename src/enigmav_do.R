enigmav_do <- function(input, ESM) {
  # The first input variable should be either "cy", meaning covariate-correction yes, or "cn", meaning covariate-correction no
  # The second input variable ESM should be either "VR" (variability ratio is computed) or "CVR" (coefficient of variation ratio is computed)
  # enigmav_do.R
  #
  # created on Fri Nov  4 17:38:56 2022
  # Philipp Homan, <philipp dot homan at bli dot uzh dot ch>
  #
  # Collaborator: Wolfgang Omlor, Finn Rabe
  #-----------------------------------------------------------------------
  #
  # load data and create timestamp for makefile

  # Check if the input is either 'cy' or 'cn'
  if (!(input %in% c("cy", "cn"))) {
    stop("Invalid input. Input should be either 'cy' or 'cn'.")
  }
  source("enigmav_load.R")
  file.create("../output/R/enigmav_do.Rout")

  # Loop through cortical, subcortical and wm measures
  # mtype <- c("cortical", "subcortical", "fa")
  mtype <- c("subcortical")

  mtype_dict <- list()
  mtype_dict$cortical <- c("CurvInd", "FoldInd", "GausCurv", "MeanCurv", "SurfAvg", "ThickAvg")
  mtype_dict$subcortical <- c("LandRvolumes")
  mtype_dict$fa <- c("FA")

  if (ESM == "VR") {
    # Load data for PBSI analysis
    iDat <- read.csv("../data/enigmav_AllIndividualDataMerged.csv")
    # Analyze and plot PBSI
    rCT <- Master_Compute_PBSI_iRS("ThickAvg", iDat)
    PBSIcA_CT <- rCT[[1]]
    PBSIpA_CT <- rCT[[2]]

    rSA <- Master_Compute_PBSI_iRS("SurfAvg", iDat)
    PBSIcA_SA <- rSA[[1]]
    PBSIpA_SA <- rSA[[2]]

    rFI <- Master_Compute_PBSI_iRS("FoldInd", iDat)
    PBSIcA_FI <- rFI[[1]]
    PBSIpA_FI <- rFI[[2]]

    rSV <- Master_Compute_PBSI_iRS("LandRvolumes", iDat)
    PBSIcA_SV <- rSV[[1]]
    PBSIpA_SV <- rSV[[2]]

    p <- plotPBSI_iBP(PBSIcA_CT, PBSIpA_CT, PBSIcA_SA, PBSIpA_SA, PBSIcA_FI, PBSIpA_FI, PBSIcA_SV, PBSIpA_SV)
    print(p)
    ggsave(filename = "../output/figures/PBSI_AllPlots.pdf", plot = p, width = 6, height = 10)

    # Analyze and plot Mean-SD relationship
    enigmav_MSDRel()
  }

  allvri_list <<- list()

  for (m in mtype) {
    print(m)
    print(mtype_dict[m])

    # At this point, main data is availabe as evp (a tibble)
    # Restrict to cortical measures for now
    corticals <- unlist(mtype_dict[m])
    evp_filt <- evp %>% dplyr::filter(measure %in% corticals)

    # This is necessary currently since there are issues with FOR and UPENN datasets
    # Simon needs to check on this
    evp_filt <- evp_filt %>% dplyr::filter(!(site %in% c("FOR2107MR", "FOR2107MS", "UPENN")))

    # Clean ROI names
    evp_filt$roic <- str_replace(evp_filt$roi, paste0("_", tolower(evp_filt$measure)), "")

    ### Take out ICV
    evp_filt <- evp_filt[evp_filt$roi != "ICV", ]
    ###

    # Get original freesurfer DK names
    dk_path <- file.path("../data", paste("dk_", m, ".txt", sep = ""))
    dk <- read.table(dk_path, col.names = "dk")
    nonrois <- c("LThickness", "RThickness", "ICV_thick")

    # Now compute summary VRs with CIs for every measure and parcel
    vri <- lapply(
      unique(evp_filt$roi),
      function(x) {
        t <- tibble(
          es = numeric(),
          ci.lb = numeric(),
          ci.ub = numeric(),
          n = numeric(),
          zval = numeric(),
          pval = numeric(),
          measure = character(),
          roi = character()
        )


        if (ESM == "VR") {
          measure <- "VR"
        } else if (ESM == "CVR") {
          measure <- "CVR"
        } else {
          stop("Invalid ESM value")
        }

        # Prepare data for random effects model
        rdat <- metafor::escalc(
          measure = measure,
          m1i = means_p, n1i = n_clean_p, sd1i = sds_p,
          m2i = means_c, n2i = n_clean_c, sd2i = sds_c,
          # data = evp_filt %>% dplyr::filter(roi == "L_bankssts_curvind"))
          data = evp_filt %>% dplyr::filter(roi == x)
        )

        ###
        mean_Age <- extend_rdat_iZ(rdat, "Age")
        mean_Sex <- extend_rdat_iZ(rdat, "Sex")
        mods_matrix <- cbind(mean_Age, mean_Sex)
        ###

        # Run random effects model to compute summary VR of
        # escalc-ed rdat
        cov_filt <- cov %>% dplyr::filter(Study %in% rdat$site)

        ###
        modX <- metafor::rma(
          yi = yi, vi = vi, data = rdat,
          method = "REML",
          slab = paste(rdat$measure),
          mods = mods_matrix,
          weighted = TRUE
        )
        ###

        mod <- metafor::rma(
          yi = yi, vi = vi, data = rdat,
          method = "REML",
          slab = paste(rdat$measure),
          weighted = TRUE
        )

        ###
        mean_age <- mean(mean_Age, na.rm = TRUE)
        mean_sex <- mean(mean_Sex, na.rm = TRUE)
        newmods <- matrix(c(mean_age, mean_sex), nrow = 1)
        # Back-transform the effect estimate and the CI
        if (input == "cy") {
          pmod <- predict(modX, newmods = newmods, trans = exp)
        } else if (input == "cn") {
          pmod <- predict(mod, trans = exp)
        }
        ###

        # Fill the tibble with estimated values
        # Also, clean ROI names (remove measure from string)
        t <- cbind(
          es = pmod$pred,
          ci.lb = pmod$ci.lb, ci.ub = pmod$ci.ub,
          n = sum(rdat$n_clean_c) + sum(rdat$n_clean_p),
          zval = mod$zval,
          pval = mod$pval,
          data.frame(measure = unique(evp_filt$measure[evp_filt$roi == x])),
          roi = str_replace(x, paste0("_", tolower(unique(evp_filt$measure[evp_filt$roi == x]))), "")
        )
      }
    )

    # Bind list to data frame and sort by effect sizes and create forest plots
    # What we have now is a tibble with summary VRs for each region; now
    # we are going to plot them per measure
    allvri <- as_tibble(do.call("rbind", vri))
    ms <- unique(allvri$measure)
    inds <- 1:length(unique(allvri$measure))

    ps <- lapply(
      inds,
      function(x) {
        dat <- filter(allvri, measure == ms[x])
        dat <- dat %>%
          arrange(desc(es))
        dat <- dat %>%
          mutate(roi = factor(roi, levels = dat$roi))
        ymin <- round(quantile(dat$ci.lb, 0.01), 2)
        ymax <- round(quantile(dat$ci.ub, 0.9), 2)
        dat$ci.lb.a <- ifelse(dat$ci.lb <= ymin, ymin, dat$ci.lb)
        dat$ci.ub.a <- ifelse(dat$ci.ub >= ymax, ymax, dat$ci.ub)
        dat$lb.arrow <- ifelse(dat$ci.lb <= ymin, 0.01, 0)
        dat$ub.arrow <- ifelse(dat$ci.ub >= ymax, 0.01, 0)
        if (ESM == "VR") {
          ylab <- paste0(ms[x], " Variability Ratio (VR)")
        } else if (ESM == "CVR") {
          ylab <- paste0(ms[x], " CVR")
        }
        fp2 <- gg_forest(dat = dat, ymin = ymin * 0.75, ymax = ymax * 1.2, ylab = ylab)[[1]]

        # adjust different plotting parameters
        if (ms[x] == "FA") {
          plot_height <- 11
          plot_width <- 7
          ggsave(
            filename = file.path("../output/", "figures", paste("enigmav_", m, "_fig.pdf", sep = "")), fp2,
            height = plot_height, width = plot_width
          )
        } else if (ms[x] == "LandRvolumes") {
          plot_height <- 5
          plot_width <- 9
          ggsave(
            filename = file.path("../output/", "figures", paste("enigmav_", m, "_fig.pdf", sep = "")), fp2,
            height = plot_height, width = plot_width
          )
        } else {
          plot_height <- 11
          plot_width <- 7
          ggsave(filename = paste0(
            "../output/figures/enigmav_cortical_figS0",
            letters[x], ".pdf"
          ), fp2, height = plot_height, width = plot_width)
        }

        if (m == "cortical") {
          vals <- dat %>%
            dplyr::filter(!roi %in% nonrois) %>%
            mutate(roi = as.character(roi))
        } else {
          vals <- dat %>%
            mutate(roi = as.character(roi))
        }
        vals <- vals[match(dk$dk, vals$roi), ]
        fn <- paste0("../output/tables/enigmav_", tolower(ms[x]), ".csv")
        df_zval <- data.frame(
          roi = vals$roi,
          zval = vals$zval
        )
        # write.table(vals$zval, fn, col.names = FALSE, row.names = FALSE)
        write.csv(df_zval, fn, col.names = TRUE, row.names = FALSE)
        rm(fn)
        p <- fp2
      }
    )

    # Connect wm labels to zval for 3d wm plot
    # system("python3 enigmav_lbl2hex.py", intern = TRUE)
    source_python("./enigmav_lbl2hex.py")
    source("enigmav_wmtracts.R")
    # Run python script that generates different plots
    # system("python3 enigmav_brains.py", intern = TRUE)
    # source_python("./enigmav_brains.py")

    # Loop through created plots and combine both forest and 3d plots in one figure
    for (i in mtype_dict[m]) {
      if (m == "cortical") {
        plot_height <- 11
        plot_width <- 7
        rel_h <- c(0.25, 1)
        imgs <- paste0("../output/figures/fsaverage_", i, "_combined.png")
      } else if (m == "subcortical") {
        plot_height <- 5.5
        plot_width <- 9
        rel_h <- c(0.5, 1)
        imgs <- paste0("../output/figures/enigmav_subcort_plot_manual.png")
      } else {
        plot_height <- 11
        plot_height <- 11
        plot_width <- 7
        rel_h <- c(0.25, 1)
        imgs <- paste0("../output/figures/enigmav_wm_fa_plot.png")
      }
      imglist <- lapply(imgs, function(x) readPNG(x))
      rglist <- lapply(imglist, function(x) rasterGrob(x, interpolate = TRUE))
      inds1 <- 1:length(rglist)
      ps2 <- lapply(inds, function(x) {
        p <- plot_grid(rglist[[x]], ps[[x]], rel_heights = rel_h, nrow = 2)

        if (ESM == "VR") {
          ggsave(paste0("../output/figures/enigmav_", m, x, "_combined.pdf"), p, height = plot_height, width = plot_width)
        } else if (ESM == "CVR") {
          ggsave(paste0("../output/figures/enigmav_", m, x, "_combined_CVR.pdf"), p, height = plot_height, width = plot_width)
        }

        pout <- p
      })
    }
    allvri_list[[m]] <<- allvri
  }
}
