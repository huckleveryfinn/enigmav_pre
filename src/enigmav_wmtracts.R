enigmav_wmtracts <- function() {

# Enable this universe
options(repos = c(
  ggseg = "https://ggseg.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

# Install some packages
# install.packages('ggsegJHU')
# install.packages('ggseg3d')
# install.packages('dplyr')
# install.packages('tidyr')
# devtools::install_github("LCBC-UiO/ggsegJHU")

# library(ggsegJHU)
# library(tidyr)
# library(ggseg3d)
# library(dplyr)

read_csv <- read.csv("../output/tables/enigmav_zval_wm_20_3d.csv")
jhu_3d$ggseg_3d[[1]]$colour <- read_csv$hex_code

zval_data <- tibble(label = read_csv$label)
zval_data <- tibble(label = jhu_3d$ggseg_3d[[1]]$region)
zval_data[, "zval"] <- read_csv$zval

# somData_aseg <- aseg_3d %>%
#   unnest(cols = ggseg_3d) %>%
#   select(label) %>%
#   filter(!grepl("Ventricle|Putamen|Amygdala", label)) %>%
#   mutate(p = seq(1, nrow(.)))

hemis <- c("left", "right")
for (h in 1:length(hemis)) {
  hemi <- hemis[h]
  fig <- ggseg3d(atlas = jhu_3d, output_dir = "./") %>%
    add_glassbrain(hemi) %>%
    pan_camera(paste(hemi, "lateral", sep = " ")) %>%
    remove_axes()
  # ggseg3d(.data = zval_data, atlas = jhu_3d, colour = "zval", text = "zval") %>% add_glassbrain("left") %>% pan_camera("right lateral") %>% remove_axes()
  # ggseg3d(.data = read_csv$zval, atlas = jhu_3d) %>% add_glassbrain("left") %>% pan_camera("right lateral") %>% remove_axes()

  # install.packages("reticulate")
  # install.packages("findpython")
  # library(findpython)
  # library(reticulate)
  # py <- find_python_cmd(minimum_version = '2.5')
  # reticulate::use_python(py, required = TRUE)
  # reticulate::py_run_string("import sys")
  # plotly::save_image(fig, width = 1400, height = 1000, file.path("output","figures",paste("enigmav_wm_",hemi,".png", sep="")))

  prt <- sample(seq(5000, 7000), 1, replace = FALSE)
  rd <- RSelenium::rsDriver(browser = "firefox", chromever = NULL, port = prt)
  #rd <- rsDriver(browser = "chrome", port = prt, verbose = T, chromever = NULL)

  # Save file and copy it to the correct location
  fn <- paste0("enigmav_wm_", hemi, ".png")
  fnpth <- file.path("../output", "figures", fn)
  plotly::export(p = fig, file = paste0("enigmav_wm_", hemi, ".png"), rd)
  print(fnpth)
  file.copy(paste0("~/Downloads/", fn), fnpth)
}
}
