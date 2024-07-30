#
# enigmav_clean.R
#
# created on Fri Nov  4 17:38:56 2022
# Philipp Homan, <philipp dot homan at bli dot uzh dot ch>
#-----------------------------------------------------------------------
source("enigmav_func.R")

# Delete final merged file (if it exists) to avoid appending of data
final_merged_file <- "../data/enigmav_AllIndividualDataMerged.csv"
if (file.exists(final_merged_file)) {
  file.remove(final_merged_file)
}

# Load aggregated data by Simon
ev <- read_csv("../preproc/all_sites/combined_group_output/table_VR_grouplevel_combined.csv")
length(unique(ev$site))

# For now get rid of any escalc measures and restrict to mean, sd, n
# Solve issue with duplicate roiname "ICV" which is associated with thickness and volume
evp <- ev %>%
  dplyr::select(
    site, modality, measure, roi, n_clean_p, means_p, sds_p,
    n_clean_c, means_c, sds_c
  )
evp$roinew <- evp$roi
evp$roinew[evp$roi == "ICV" & evp$measure == "ThickAvg"] <- "ICV_thick"
evp$roi <- evp$roinew
evp <- evp %>% dplyr::select(-roinew)

# Save data file to disk
write_csv(evp, file = "../data/enigmav_data.csv")

# Visualize
## p1 <- ggplot(data=evp %>% filter(measure=="ThickAvg",
##                                  !roi=="ICV_thick",
##                                  !site %in% c("FOR2107MR", "FOR2107MS")),
##              aes(x=means_p)) +
##   facet_wrap(~roi) +
##   geom_histogram()

## p2 <- ggplot(data=evp %>% filter(measure=="SurfAvg",
##                                  !roi=="LSurfArea",
##                                  !roi=="RSurfArea",
##                                  !site %in% c("FOR2107MR", "FOR2107MS")),
##              aes(x=means_p)) +
##   facet_wrap(~roi) +
##   geom_histogram()

## p3 <- ggplot(data=evp %>% filter(measure=="LandRvolumes",
##                                  !roi=="ICV",
##                                  !site %in% c("FOR2107MR", "FOR2107MS")),
##              aes(x=means_p)) +
##   facet_wrap(~roi) +
##   geom_histogram()


## # Find extreme values
## # as.data.frame(evp %>% filter(measure=="ThickAvg", means_p > 100))
## # as.data.frame(evp %>% filter(measure=="SurfAvg", means_p > 10000))

# Deal with single sites and subject level data
# Note thatpreproc/all_sites has
# been cleaned so that it only includes site subfolders

# First, obtain a vector with all sites
sites <- list.dirs(
  path = "../preproc/all_sites/",
  full.names = FALSE,
  recursive = FALSE
)

measures <- c(
  "CorticalMeasuresENIGMA_CurvInd",
  "CorticalMeasuresENIGMA_FoldInd",
  "CorticalMeasuresENIGMA_GausCurv",
  "CorticalMeasuresENIGMA_MeanCurv",
  "CorticalMeasuresENIGMA_SurfAvg",
  "CorticalMeasuresENIGMA_ThickAvg",
  "SubcorticalMeasuresENIGMA_LandRvolumes"
  # "DTI_FA_allROIs_grouplevel"
)

measures_short <- unlist(str_split(measures, "_"))[c(2, 4, 6, 8, 10, 12, 14)]

# Each data file needs to be transformed to long format
prepth <- "../preproc/all_sites/"

dlist <- vector("list", length = length(sites) * length(measures))
VecX <- list()
count <- 0
count2 <- 0
countY <- 0
n_removed_df <- data.frame(file=character(), n_removed=integer(), stringsAsFactors=FALSE)

all_dataX <- list()

for (i in 1:length(sites)) {
  
  print(paste0("Processing site ", i, " out of ", length(sites), ": ", sites[i]))
  
  dlj1f <- paste0(prepth, sites[i], "/", sites[i], "_allmerged.csv")
  
  # At the start of processing each site, delete existing CSV file
  if (file.exists(dlj1f)) {
    file.remove(dlj1f)
  }
  
  cf <- paste0(prepth, sites[i], "/", "Covariates.csv")
  if (file.exists(cf)) {
    covars <- suppressMessages(read_csv(cf, show_col_types = FALSE))
    countY = countY + 1
    if (sites[i] == "UNINA"){
      # Identify problematic rows with "Clotiapine", "Amisulpride", or "Amisilpride"
      idx <- which(grepl("Clotiapine|Amisulpride|Amisilpride", covars$"ANTIPSYCHOTIC_TYPE_(typical/atypical/both)"))
      if(length(idx) > 0){
        # Create a new variable to hold the values from columns 7 to 14 of covars
        new_data <- covars[idx, 7:14]
        # Change the data types of the new_data to match with one column before in covars
        for(col in seq_along(new_data)) {
          new_data[, col] <- get(paste0("as.", typeof(covars[idx, col + 5])))(new_data[, col])
        }
        # Ensure that column names of new_data and corresponding columns in covars are the same
        colnames(new_data) <- colnames(covars[, 7:14])
        # Create a new dataframe by removing problematic rows from covars
        covars_new <- covars[-idx, ]
        # Add new_data to corresponding rows in the new dataframe
        for(zzz in seq_along(idx)) {
          new_row <- data.frame(covars[idx[zzz], 1:5], new_data[zzz, ], covars[idx[zzz], 14:ncol(covars)])
          colnames(new_row) <- colnames(covars) # ensure column names match
          covars_new <- rbind(covars_new, new_row)
        }
        covars_new <- covars_new[order(as.numeric(rownames(covars_new))), ]
      }
      covars_new <- covars_new[, -14]
      # Modify "AP" column
      covars_new$AP <- ifelse(covars_new$Dx == 0, 0,
                              ifelse(covars_new$Dx == 1 & covars_new$`ANTIPSYCHOTIC_TYPE_(typical/atypical/both)` == "atypical", 2,
                                     ifelse(covars_new$Dx == 1 & covars_new$`ANTIPSYCHOTIC_TYPE_(typical/atypical/both)` == "typical", 3,
                                            ifelse(covars_new$Dx == 1 & covars_new$`ANTIPSYCHOTIC_TYPE_(typical/atypical/both)` == "both", 4, NA))))
      
      if(!is.numeric(covars_new$AP)) covars_new$AP <- as.numeric(covars_new$AP)
      covars_new$HAND <- ifelse(covars_new$HAND == "dx", 0,
                                ifelse(covars_new$HAND == "sx", 1, NA))
      if(!is.numeric(covars_new$HAND)) covars_new$HAND <- as.numeric(covars_new$HAND)
      # Convert "CPZ", "AO" and "DURILL" to numeric
      numeric_cols <- c("CPZ", "DURILL", "AO")
      for (col in numeric_cols) {
        if (col %in% names(covars_new)) {
          covars_new[[col]] <- ifelse(covars_new[[col]] == "", NA, covars_new[[col]])
          covars_new[[col]] <- as.numeric(as.character(covars_new[[col]]))
        }
      }     
      covars <- covars_new
    }
    # Transform Sex column
    covars <- covars %>%
      dplyr::mutate(Sex = ifelse(Sex == "M", 1,
                                 ifelse(Sex == "F", 2, Sex)))
    
    n_before <- nrow(covars)
    covars <- covars %>%
      distinct(SubjID, .keep_all = TRUE)
    n_after <- nrow(covars)
    n_removed <- n_before - n_after
    if(n_removed > 0){n_removed_df <- rbind(n_removed_df, data.frame("file"=sites[i], "n_removed"=n_removed, stringsAsFactors=FALSE))}
  }
  
  allowed_cols <- c('SubjID', 'Dx', 'Age', 'Sex', 'AP', 'CPZ', 'AO', 'DURILL', 'PANSSTOT', 'PANSSPOS', 'PANSSNEG', 'HAND')
  if(exists("covars")){covars <- covars %>% 
    dplyr::select(which(colnames(covars) %in% allowed_cols))
  if("Age" %in% names(covars)){
    covars <- covars %>%
      filter((Age >= 18 & Age <= 65))
  }
  if("AO" %in% names(covars)){
    covars <- covars %>%
      filter((AO >= 18 & AO <= 65) | is.na(AO))
  }
  }
  countX <- 0
  for (j in 1:length(measures)) {
    f <- paste0(prepth, sites[i], "/", measures[j], ".csv")
    if (file.exists(f)) {
      # Read file and transform from wide to long
      count <- count + 1
      countX <- countX + 1
      VecX[[countY]] <- countX
      d <- suppressMessages(read_csv(f, show_col_types = FALSE))
      d1 <- d %>%
        dplyr::mutate(newid = seq.int(nrow(d)))
      dl <- d1 %>%
        gather(key = roi, value = meanval, c(-1, -length(d1))) %>%
        mutate(
          measure = measures_short[j],
          site = sites[i],
          id = paste0(sites[i], "_", as.character(newid))
        )
      # Join with covariates if possible
      if (exists("covars")) {
        covars$Sex <- as.character(covars$Sex)
        dlj <- dl %>% left_join(covars, by = "SubjID")
      }
      dlj1 <- dlj %>%
        dplyr::select(-SubjID)
      
      dlj1 <- dlj1 %>%
        dplyr::mutate(Sex = ifelse(Sex == "M", 1,
                                   ifelse(Sex == "F", 2, Sex)))
      dlist[[count]] <- dlj1
      dlj1f <- paste0(prepth, sites[i], "/", sites[i], "_allmerged.csv")

      if (file.exists(dlj1f)) {
        existing_data <- suppressMessages(read_csv(dlj1f))
        existing_data$Sex <- as.character(existing_data$Sex)
        # Convert the specified columns to numeric
        numeric_cols <- c("CPZ", "DURILL", "AP", "AO")
        for (col in numeric_cols) {
          if (col %in% names(existing_data)) {
            existing_data[[col]] <- ifelse(existing_data[[col]] == "", NA, existing_data[[col]])
            existing_data[[col]] <- as.numeric(as.character(existing_data[[col]]))
          }
          if (col %in% names(dlj1)) {
            dlj1[[col]] <- ifelse(dlj1[[col]] == "", NA, dlj1[[col]])
            dlj1[[col]] <- as.numeric(as.character(dlj1[[col]]))
          }
        } 
        new_data <- suppressMessages(bind_rows(existing_data, dlj1))
      } else {
        new_data <- dlj1
      }
      print(paste0("Appending data frame ", count, " to list"))
      all_dataX[[count]] <- new_data
      print(paste0("Written file: ", dlj1f))
      write_csv(new_data, dlj1f)  
    }
  }
  processed_covars_file <- paste0(prepth, sites[i], "/Covariates_processed.csv")
  if (file.exists(processed_covars_file)) {
    file.remove(processed_covars_file)
  }
  if(exists("covars")){write_csv(covars, file = paste0(prepth, sites[i], "/Covariates_processed.csv"))
  rm(covars)}
  }
print("Binding all data frames in list together")
VecX <- cumsum(unlist(VecX))
#all_data_df <- bind_rows(all_dataX)
all_data_df <- bind_rows(all_dataX[VecX])
write_csv(all_data_df, file = "../data/enigmav_AllIndividualDataMerged.csv")

CountMeasures_iZ()
Meta_ScanCovariates_iZ()
CountRows4Sites_iZ()
print("Number of removed duplicate IDs:")
print(n_removed_df)