
# CODE FOR ORGANIZING SRS AND SLEEP DATA FROM RAW, INCLUDING IMPUTING MISSING DATA AND RESIDUALIZING 2/27/2025 -- Updated 4/2/2025 with new data drop

rm(list = ls())
library(tidyverse)
library(mice)
library(OpenMx)

# Read in raw data
vcu_07032024_raw <- read_csv("./002_ASDdata07032024.csv")
vcu_03192024_raw <- read_csv("./001_ASDdata03192024.csv")
vcu_09182024_raw <- read_csv("./003_ASDdata09182024.csv")
vcu_11062024_raw <- read_csv("./004_ASDdata11062024.csv")
vcu_01222025_raw <- read_csv("./005_ASDdata01222025.csv")
vcu_04022025_raw <- read_csv("./006ASDdata04022025.csv")

# Remove extra column
vcu_07032024_raw <- vcu_07032024_raw %>%
  select(-notes)
vcu_09182024_raw <- vcu_09182024_raw %>%
  select(-notes)
vcu_11062024_raw <- vcu_11062024_raw %>%
  select(-notes)
vcu_01222025_raw <- vcu_01222025_raw %>%
  select(-notes)
vcu_04022025_raw <- vcu_04022025_raw %>%
  select(-notes)


# Merging data frames together
vcu_merged_raw <- rbind(vcu_07032024_raw, vcu_03192024_raw, vcu_09182024_raw, vcu_11062024_raw, vcu_01222025_raw, vcu_04022025_raw)

# Choosing columns of interest
vcu_dem_srs_sds_raw <- vcu_merged_raw %>%
  select(efamno, pair, eid, zygo, msex, Twins_race, dem_age_twn_yrs, dem_twina_sex, dem_twinb_sex, hlth_twina_condition___3, hlth_twina_condition___13, hlth_twina_conditiono, hlth_twinb_condition___3, hlth_twinb_condition___13, hlth_twinb_conditiono, srs_twa_1:srs_twb_65, sds_twa_duration:sds_twb24)

vcu_dem_srs_sds_raw <- vcu_dem_srs_sds_raw %>%
  rename(FID = efamno) %>% # family ID
  rename(PNUM = pair) %>% # pair number in family
  rename(TID = eid) %>% # twin pair ID
  rename(race = Twins_race)


# ------------------------------------------------------------------------------------------------------------------
# RECODE RACE VARIABLE to match the five category system in ABCD for ease of analysis
# ------------------------------------------------------------------------------------------------------------------

# White = 0, African American = 1, Hispanic = 2, Asian or Pacific Islander = 3, Other/Mixed Race = 4, Unknown/NA = 5
vcu_dem_srs_sds <- vcu_dem_srs_sds_raw %>%
  mutate(race = case_when(
    race == "White" ~ 0,
    race == "African American" ~ 1,
    race == "Hispanic" ~ 2,
    race == "Asian or Pacific Islander" ~ 3,
    race == "Unknown" ~ 5,
    race == "Other" ~ 4,
    race == "Mixed Race" ~ 4,
  ))

vcu_dem_srs_sds$race[is.na(vcu_dem_srs_sds$race)] <- 5

# Recode race/ethnicity variables
vcu_dem_srs_sds$isCA <- 0
vcu_dem_srs_sds$isCA[vcu_dem_srs_sds$race == 0]     <- 1 # White/caucasian
vcu_dem_srs_sds$isAA <- 0
vcu_dem_srs_sds$isAA[vcu_dem_srs_sds$race == 1]     <- 1 # AfAm/Black
vcu_dem_srs_sds$isHI <- 0
vcu_dem_srs_sds$isHI[vcu_dem_srs_sds$race == 2]  <- 1 # Hispanic
vcu_dem_srs_sds$isAS <- 0
vcu_dem_srs_sds$isAS[vcu_dem_srs_sds$race == 3]     <- 1 # Asian or Pacific Islander
vcu_dem_srs_sds$isMI <- 0
vcu_dem_srs_sds$isMI[vcu_dem_srs_sds$race == 4]     <- 1 # Mixed or other race
vcu_dem_srs_sds$isOT <- 0
vcu_dem_srs_sds$isOT[vcu_dem_srs_sds$race == 5]     <- 1 # Other/unknown

# Count the number of 1s in each race category
race_counts <- colSums(vcu_dem_srs_sds[, c("isCA", "isAA", "isHI", "isAS", "isMI", "isOT")])
print(race_counts)

# Calculating reverse scores for twin a
vcu_dem_srs_sds_scored <- vcu_dem_srs_sds %>%
  mutate(srs_twa_3 = 5 - srs_twa_3) %>%
  mutate(srs_twa_7 = 5 - srs_twa_7) %>%
  mutate(srs_twa_11 = 5 - srs_twa_11) %>%
  mutate(srs_twa_12 = 5 - srs_twa_12) %>%
  mutate(srs_twa_15 = 5 - srs_twa_15) %>%
  mutate(srs_twa_17 = 5 - srs_twa_17) %>%
  mutate(srs_twa_21 = 5 - srs_twa_21) %>%
  mutate(srs_twa_22 = 5 - srs_twa_22) %>%
  mutate(srs_twa_26 = 5 - srs_twa_26) %>%
  mutate(srs_twa_32 = 5 - srs_twa_32) %>%
  mutate(srs_twa_38 = 5 - srs_twa_38) %>%
  mutate(srs_twa_40 = 5 - srs_twa_40) %>%
  mutate(srs_twa_43 = 5 - srs_twa_43) %>%
  mutate(srs_twa_45 = 5 - srs_twa_45) %>%
  mutate(srs_twa_48 = 5 - srs_twa_48) %>%
  mutate(srs_twa_52 = 5 - srs_twa_52) %>%
  mutate(srs_twa_55 = 5 - srs_twa_55)

# Calculating reverse scores for twin b
vcu_dem_srs_sds_scored <- vcu_dem_srs_sds_scored %>%
  mutate(srs_twb_3 = 5 - srs_twb_3) %>%
  mutate(srs_twb_7 = 5 - srs_twb_7) %>%
  mutate(srs_twb_11 = 5 - srs_twb_11) %>%
  mutate(srs_twb_12 = 5 - srs_twb_12) %>%
  mutate(srs_twb_15 = 5 - srs_twb_15) %>%
  mutate(srs_twb_17 = 5 - srs_twb_17) %>%
  mutate(srs_twb_21 = 5 - srs_twb_21) %>%
  mutate(srs_twb_22 = 5 - srs_twb_22) %>%
  mutate(srs_twb_26 = 5 - srs_twb_26) %>%
  mutate(srs_twb_32 = 5 - srs_twb_32) %>%
  mutate(srs_twb_38 = 5 - srs_twb_38) %>%
  mutate(srs_twb_40 = 5 - srs_twb_40) %>%
  mutate(srs_twb_43 = 5 - srs_twb_43) %>%
  mutate(srs_twb_45 = 5 - srs_twb_45) %>%
  mutate(srs_twb_48 = 5 - srs_twb_48) %>%
  mutate(srs_twb_52 = 5 - srs_twb_52) %>%
  mutate(srs_twb_55 = 5 - srs_twb_55)

vcu_dem_srs_sds_scored <- vcu_dem_srs_sds_scored %>%
  mutate(
    zygo_group = case_when(
      zygo == "MZ" & dem_twina_sex == 1 & dem_twinb_sex == 1 ~ 1,  # Male MZ
      zygo == "MZ" & dem_twina_sex == 2 & dem_twinb_sex == 2 ~ 2,  # Female MZ
      zygo == "DZ" & dem_twina_sex == 1 & dem_twinb_sex == 1 ~ 3,  # Male DZ
      zygo == "DZ" & dem_twina_sex == 2 & dem_twinb_sex == 2 ~ 4,  # Female DZ
      zygo == "DZ" & dem_twina_sex != dem_twinb_sex ~ 5          # Opposite-sex DZ
    )
  )

# Making separate twin a and twin b data frames
vcu_dem_srs_sds_scored_twa <- vcu_dem_srs_sds_scored %>%
  select(FID, PNUM, TID, zygo, race, dem_age_twn_yrs, dem_twina_sex, hlth_twina_condition___3, hlth_twina_condition___13, isCA:isOT, srs_twa_1:srs_twa_65, sds_twa_duration:sds_twa24, zygo_group)

vcu_dem_srs_sds_scored_twb <- vcu_dem_srs_sds_scored %>%
  select(FID, PNUM, TID, zygo, race, dem_age_twn_yrs, dem_twinb_sex, hlth_twinb_condition___3, hlth_twinb_condition___13, isCA:isOT, srs_twb_1:srs_twb_65, sds_twb_duration:sds_twb24, zygo_group) %>%
  rename(sds_twb_sleeponset = sds_twa_sleeponset_2)

# Rename columns
library(stringr)
vcu_dem_srs_sds_scored_twa <- vcu_dem_srs_sds_scored_twa %>%
  rename_with(~ str_replace(., "^srs_twa_", "srs_"), starts_with("srs_twa_")) %>%
  rename_with(~ str_replace(., "^sds_twa_", "sds_"), starts_with("sds_twa_")) %>%
  rename_with(~ str_replace(., "^sds_twa", "sds_"), starts_with("sds_twa")) %>%
  rename(age = dem_age_twn_yrs,
         sex = dem_twina_sex,
         dx_asd = hlth_twina_condition___3,
         dx_aspergers =  hlth_twina_condition___13) %>%
  mutate(orig_id = 8) # originally twin a
vcu_dem_srs_sds_scored_twb <- vcu_dem_srs_sds_scored_twb %>%
  rename_with(~ str_replace(., "^srs_twb_", "srs_"), starts_with("srs_twb_")) %>%
  rename_with(~ str_replace(., "^sds_twb_", "sds_"), starts_with("sds_twb_")) %>%
  rename_with(~ str_replace(., "^sds_twb", "sds_"), starts_with("sds_twb")) %>%
  rename(age = dem_age_twn_yrs,
         sex = dem_twinb_sex,
         dx_asd = hlth_twinb_condition___3,
         dx_aspergers =  hlth_twinb_condition___13) %>%
  mutate(orig_id = 9) # originally twin b

# Combining df to get all of the twins in one long df (each twin has unique row)
vcu_long <- rbind(vcu_dem_srs_sds_scored_twa, vcu_dem_srs_sds_scored_twb) %>%
  mutate(IID = row_number())

# write.csv(vcu_long, "./Amanda_analysis/01_tidy_data/vcu_long.csv", row.names = FALSE)

vcu_long <- vcu_long %>%
  mutate(srs_tot_raw = rowSums(select(., srs_1, srs_2, srs_3, srs_4, srs_5, srs_6, srs_7, srs_8, srs_9, srs_10, srs_11, srs_12, srs_13, srs_14, srs_15, srs_16, srs_17, srs_18, srs_19, srs_20, srs_21, srs_22, srs_23, srs_24, srs_25, srs_26, srs_27, srs_28, srs_29, srs_30, srs_31, srs_32, srs_33, srs_34, srs_35, srs_36, srs_37, srs_38, srs_39, srs_40, srs_41, srs_42, srs_43, srs_44, srs_45, srs_46, srs_47, srs_48, srs_49, srs_50, srs_51, srs_52, srs_53, srs_54, srs_55, srs_56, srs_57, srs_58, srs_59, srs_60, srs_61, srs_62, srs_63, srs_64, srs_65)))

vcu_long <- vcu_long %>%
  mutate(srs_stot_raw = rowSums(select(., srs_6, srs_15, srs_16, srs_18, srs_24, srs_29, srs_35, srs_37, srs_39, srs_42, srs_58)))

vcu_long <- vcu_long %>%
  mutate(sds_tot_raw = rowSums(select(., sds_duration:sds_24)))

# Make df with only srs and sds
vcu_items <- vcu_long %>%
  select(IID, srs_1:srs_65, sds_duration:sds_24)

# Function to calculate the percentage of missing responses per participant
pMiss <- function(x) { sum(is.na(x)) / length(x) * 100 }

# Apply the function row-wise to compute missing percentage for each participant
vcu_items$missing_percent <- apply(vcu_items, 1, pMiss)

# Remove participants that are missing 5% or more of data
vcu_cleaned <- vcu_items %>% filter(missing_percent < 5)

# Select only the IDs from vcu_cleaned
ids_to_keep <- vcu_cleaned %>% select(IID)

# Filter vcu_wide_clean to keep only individuals with minimal missing data
vcu_long_clean <- vcu_long %>% semi_join(ids_to_keep, by = "IID")

# Separate data for srs_tot and srs_stot
vcu_srs <- vcu_long_clean %>%
  select(IID, FID, PNUM, TID, zygo, race, age, sex, dx_asd, dx_aspergers, isCA:isOT, srs_1:srs_65, orig_id)
vcu_sds <- vcu_long_clean %>%
  select(IID, FID, PNUM, TID, zygo, race, age, sex, dx_asd, dx_aspergers, isCA:isOT, sds_duration:sds_24, orig_id)

# Impute srs and sds
pred_matrix_1 <- quickpred(vcu_srs, exclude = c("FID", "PNUM", "TID"))  # Adjust as needed
imp_srs <- mice(vcu_srs, method = "cart", pred = pred_matrix_1, m = 1, maxit = 10, seed = 123)

pred_matrix_2 <- quickpred(vcu_sds, exclude = c("FID", "PNUM", "TID"))  # Adjust as needed
imp_sds <- mice(vcu_sds, method = "cart", pred = pred_matrix_2, m = 1, maxit = 10, seed = 123)

vcu_srs_imp <- complete(imp_srs)
vcu_sds_imp <- complete(imp_sds)

sum(is.na(vcu_srs))  # Count missing values before imputation
sum(is.na(complete(imp_srs)))  # Should be 0 if all missing values were imputed

sum(is.na(vcu_sds))  # Count missing values before imputation
sum(is.na(complete(imp_sds)))  # Should be 0 if all missing values were imputed

# Score forms now that missing data is imputed
vcu_srs_imp <- vcu_srs_imp %>%
  mutate(srs_tot_imp = rowSums(select(., srs_1, srs_2, srs_3, srs_4, srs_5, srs_6, srs_7, srs_8, srs_9, srs_10, srs_11, srs_12, srs_13, srs_14, srs_15, srs_16, srs_17, srs_18, srs_19, srs_20, srs_21, srs_22, srs_23, srs_24, srs_25, srs_26, srs_27, srs_28, srs_29, srs_30, srs_31, srs_32, srs_33, srs_34, srs_35, srs_36, srs_37, srs_38, srs_39, srs_40, srs_41, srs_42, srs_43, srs_44, srs_45, srs_46, srs_47, srs_48, srs_49, srs_50, srs_51, srs_52, srs_53, srs_54, srs_55, srs_56, srs_57, srs_58, srs_59, srs_60, srs_61, srs_62, srs_63, srs_64, srs_65)))

vcu_srs_imp <- vcu_srs_imp %>%
  mutate(srs_stot_imp = rowSums(select(., srs_6, srs_15, srs_16, srs_18, srs_24, srs_29, srs_35, srs_37, srs_39, srs_42, srs_58)))

vcu_sds_imp <- vcu_sds_imp %>%
  mutate(sds_tot_imp = rowSums(select(., sds_duration:sds_24))) %>%
  select(IID, sds_duration:sds_tot_imp)

# Merge with the SRS dataset
vcu_long_imp <- left_join(vcu_srs_imp, vcu_sds_imp, by = "IID") %>%
  select(IID:isOT, srs_tot_imp, srs_stot_imp, sds_tot_imp) %>%
  rename(srs_imp = srs_tot_imp) %>%
  rename(ssrs_imp = srs_stot_imp) %>%
  rename(sds_imp = sds_tot_imp)

# write.csv(vcu_long_imp, "./Amanda_analysis/01_tidy_data/vcu_long_imp.csv", row.names = FALSE)

#vcu_long_imp <- read_csv("./Amanda_analysis/01_tidy_data/vcu_long_imp.csv")

# ------------------------------------------------------------------------------------------------------------------
# RESIDUALIZE DATA
# ------------------------------------------------------------------------------------------------------------------

vcu25l_imp <- vcu_long_imp

# Identify variables of interest
vcuVars <- names(vcu25l_imp[, c(17:19)])  # Variables for regression

# Perform regression and store residuals
for (i in vcuVars) {
  # Fit linear model
  temp <- lm(vcu25l_imp[[i]] ~ age + isCA + isAA + isHI + isAS + isMI + isOT, 
             data = vcu25l_imp, 
             na.action = na.exclude)
  
  # Store residuals in a new column
  vcu25l_imp[[paste0('r', i)]] <- residuals(temp)
}

# ------------------------------------------------------------------------------------------------------------------
# REMOVE OUTLIERS
# ------------------------------------------------------------------------------------------------------------------

# Function to standardize variables
myStan <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Select the variables to standardize
rVars <- names(vcu25l_imp)[20:22] # Standardizing residualized data

# Standardize and handle outliers
for (var in rVars) {
  # Standardize the column
  standardized <- myStan(vcu25l_imp[[var]])
  
  # Replace outliers (values beyond ±4 SD) with NA
  standardized[standardized < -4 | standardized > 4] <- NA
  
  # Assign back to the corresponding column in vcu25l
  vcu25l_imp[[var]] <- standardized
}

# Count the number of NAs per column
na_counts_per_column <- colSums(is.na(vcu25l_imp))
print(na_counts_per_column)

# Remove rows with any NA values (outliers)
vcu25lrs_imp <- na.omit(vcu25l_imp) #long, residualized, standardized

# ------------------------------------------------------------------------------------------------------------------
# RANDOMIZE TWINS AND MAKE WIDE DATA FRAME (each family has unique row)
# ------------------------------------------------------------------------------------------------------------------

# Identify unmatched TIDs (those that appear only once). These are now unmatched because their twin had an outlier score
unmatched_TIDs <- vcu25lrs_imp %>%
  group_by(TID) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(TID) %>%
  distinct()

# Remove rows where TID has no pair (i.e., those in unmatched_TIDs)
vcu25lrso_imp <- vcu25lrs_imp %>% # NAs and outliers removed
  filter(!TID %in% unmatched_TIDs$TID)

# Check that every twin has a pair (the result should be empty df)
vcu25lrso_imp %>%
  group_by(TID) %>%
  summarise(n = n()) %>%
  filter(n != 2)

set.seed(123)  # Set the seed to 123
vcu25lrso_imp <- vcu25lrso_imp %>%
  group_by(TID) %>%
  mutate(ah_id_new = sample(c("1", "2"))) %>%
  ungroup()

# Check that each pair has a T1 and a T2 (the result should be empty df)
vcu25lrso_imp %>%
  group_by(TID, ah_id_new) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n != 1)


# Make separate dfs for T1 and T2
vcu25_imp_T1 <- vcu25lrso_imp %>%
  filter(ah_id_new == "1") %>%
  select(IID:zygo, sex, srs_imp:ah_id_new, -PNUM, -FID) %>%
  rename(
    IID_T1 = IID,
    sex_T1 = sex,
    srs_imp_T1 = srs_imp,
    sds_imp_T1 = sds_imp,
    ssrs_imp_T1 = ssrs_imp,
    rsrs_imp_T1 = rsrs_imp,
    rsds_imp_T1 = rsds_imp,
    rssrs_imp_T1 = rssrs_imp
  )

vcu25_imp_T2 <- vcu25lrso_imp %>%
  filter(ah_id_new == "2") %>%
  select(IID:zygo, sex, srs_imp:ah_id_new, -PNUM, -FID) %>%
  rename(
    IID_T2 = IID,
    sex_T2 = sex,
    srs_imp_T2 = srs_imp,
    sds_imp_T2 = sds_imp,
    ssrs_imp_T2 = ssrs_imp,
    rsrs_imp_T2 = rsrs_imp,
    rsds_imp_T2 = rsds_imp,
    rssrs_imp_T2 = rssrs_imp
  )

# Putting dfs together
vcu25w_imp <- merge(vcu25_imp_T1, vcu25_imp_T2, by = "TID") %>%
  rename(zygo = zygo.x) %>%
  select(-zygo.y, -ah_id_new.x, -ah_id_new.y) %>%
  mutate(
    zygo_group = case_when(
      zygo == "MZ" & sex_T1 == 1 & sex_T2 == 1 ~ 1,  # Male MZ
      zygo == "MZ" & sex_T1 == 2 & sex_T2 == 2 ~ 2,  # Female MZ
      zygo == "DZ" & sex_T1 == 1 & sex_T2 == 1 ~ 3,  # Male DZ
      zygo == "DZ" & sex_T1 == 2 & sex_T2 == 2 ~ 4,  # Female DZ
      zygo == "DZ" & sex_T1 != sex_T2 ~ 5          # Opposite-sex DZ
    )
  )

# Save df
# write.csv(vcu25w_imp, "./Amanda_analysis/01_tidy_data/vcu25w_imp.csv", row.names = FALSE)
write.csv(vcu25w_imp, "./vcu25w_imp_paper.csv", row.names = FALSE)

