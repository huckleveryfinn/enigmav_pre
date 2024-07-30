#
# enigmav_func.R
#
# created on Fri Nov  4 17:38:56 2022
# Philipp Homan, <philipp dot homan at bli dot uzh dot ch>
#-----------------------------------------------------------------------
#
# common libraries


libs <- c(
  "tidyr",
  "dplyr",
  "devtools",
  "png",
  "grid",
  "broom",
  "represearch",
  "metafor",
  "ggplot2",
  "tidyverse",
  "cowplot",
  "knitr",
  "rmarkdown",
  "papaja",
  "kableExtra",
  "ggsegJHU",
  "ggseg3d",
  "ggsegJHU",
  "readr",
  "reticulate",
  "findpython",
  "ppcor",
  "gridExtra",
  "RSelenium",
  "Hmisc"
)

if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(char = libs)

# use_python("/Users/frabe/anaconda3/bin/python")

source("Plot_CorrMSD.R")
source("PlotCorrMSD_SumStats.R")
source("Analyze_MSDRel.R")
source("enigmav_MSDRel.R")

source("Compute_PBSIag_iZ.R")
source("Compute_PBSIspec_iRS.R")
source("Meta_Compute_PBSI_iRS.R")
source("Master_Compute_PBSI_iRS.R")
source("Plot_PBSI_iZ.R")

source("enigmav_wmtracts.R")
source("CountMeasures_iZ.R")
source("ScanCovariates_iZ.R")
source("Meta_ScanCovariates_iZ.R")
source("CountRows4Sites_iZ.R")
source("extend_rdat_iZ.R")

source("plotPBSIX10.R")
source("plotPBSI_iBP.R")
source("Create_ErrorbarPlot.R")

# Forest plot with ggplot
#-----------------------------------------------------------------------
gg_forest <- function(dat, ymin, ymax, type = "VR",
                      xlab = "Region of interest",
                      ylab = "Variability Ratio (VR)",
                      cilab = "VR [95% CI]") {
  # if (type == "VR") {
  #  ylab <- "Variability Ratio (VR)"
  #  cilab <- "VR [95% CI]"
  # } else {
  #  ylab <- "Coefficient of Variation Ratio (CVR)"
  #  cilab <- "CVR [95% CI]"
  # }
  # ymax <- ymax * 0.8

  fp1 <- ggplot(dat, aes(x = roi, y = es)) +
    # ggtitle(element_text(unique(dat$measure), face="bold")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_point(size = 2) +
    geom_segment(
      aes(
        x = roi, xend = roi, y = es,
        yend = ci.ub.a
      ),
      size = 0.4,
      arrow = arrow(
        length = unit(dat$ub.arrow, "npc"),
        ends = "last", type = "closed"
      )
    ) +
    geom_segment(
      aes(
        x = roi, xend = roi, y = es,
        yend = ci.lb.a
      ),
      size = 0.4,
      arrow = arrow(
        length = unit(dat$lb.arrow, "npc"),
        ends = "last", type = "closed"
      )
    ) +
    scale_shape_manual(values = c(19, 18)) +
    scale_color_manual(values = c("black", "white")) +
    scale_size(range = c(1, 2)) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      # axis.line.x = element_line(),
      # axis.ticks.x = element_line(),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 7, color = "black")
    ) +
    ylab(ylab) +
    xlab(xlab) +
    ylim(ymin, ymax * 1.08)
  fp2 <- fp1 +
    annotate("text",
      y = ymax * 1.08, x = nrow(dat) + 1.5,
      label = paste0(cilab), hjust = 1, fontface = "bold"
    ) +
    annotate("text",
      y = ymax * 1.08, x = 1:nrow(dat), hjust = 1,
      label = paste0(
        sprintf("%.2f", round(dat$es, 2)),
        " [",
        sprintf("%.2f", round(dat$ci.lb, 2)),
        ", ",
        sprintf("%.2f", round(dat$ci.ub, 2)),
        "]"
      ), size = 2.5
    ) +
    annotate("text",
      y = 0.99, x = nrow(dat) + 1.5,
      label = "Greater in controls", hjust = 1,
      fontface = "bold"
    ) +
    annotate("text",
      y = 1.01, x = nrow(dat) + 1.5,
      label = "Greater in cases",
      hjust = 0, fontface = "bold"
    ) +
    annotate("text",
      y = ymin * 1.00, x = nrow(dat) + 1.5,
      label = "N", hjust = 1, fontface = "bold"
    ) +
    annotate("text",
      y = ymin * 1.00, x = 1:nrow(dat), label = paste0(dat$n),
      hjust = 1, size = 2.5
    )
  fp2 <- fp2 + coord_flip(clip = "off")
  return(list(fp2))
}
