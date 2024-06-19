
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(stringr)
library(data.table)

# Data --------------------------------------------------------------------

all_data <- read_csv("./labels_new_final_dated.csv")
all_data$date_of_study <- as.Date(all_data$date_of_study, format = "%m/%d/%Y")

## adni data ----------------------------------------------------------

adni <- all_data[nchar(as.character(all_data$patient_original_id)) == 10, ]
adni <- adni %>% mutate(RID = str_sub(patient_original_id, start = -4))
adni_rids <- adni$RID
adni_rids[str_sub(adni_rids, end = 1) == 0] <- str_sub(adni_rids[str_sub(adni_rids, end = 1) == 0], start = 2)
adni$RID <- adni_rids
adni_ids <- unique(adni$patient_original_id)
adni_rids <- unique(adni$RID)

ad_dxsum <- read_csv("../useful clinical data/ad_mci_hc_dxsum.csv")
ad_dxsum <- ad_dxsum %>% filter(PTID %in% adni_ids)
ad_dxsum <- ad_dxsum[order(ad_dxsum$PTID, decreasing = FALSE, na.last = FALSE),]

ad_mmse <- read_csv("../useful clinical data/ad_mci_hc_mmse_trimmed.csv")
ad_mmse <- ad_mmse %>% filter(RID %in% adni_rids)
ad_mmse$USERDATE = as.Date(ad_mmse$USERDATE, format = "%m/%d/%Y")

ad_cdr <- read_csv("../useful clinical data/ad_mci_hc_cdr_trimmed.csv")
ad_cdr <- ad_cdr %>% filter(RID %in% adni_rids)
ad_cdr$USERDATE = as.Date(ad_cdr$USERDATE, format = "%m/%d/%Y")

ad_demog <- read_csv("../useful clinical data/ad_mci_hc_demographics_trimmed.csv")
ad_demog <- ad_demog %>% filter(RID %in% adni_rids)
ad_demog$USERDATE = as.Date(ad_demog$USERDATE, format = "%m/%d/%Y")

ad_gds <- read_csv("../useful clinical data/ad_mci_hc_gds_trimmed.csv")
ad_gds <- ad_gds %>% filter(RID %in% adni_rids)
ad_gds$USERDATE = as.Date(ad_gds$USERDATE, format = "%m/%d/%Y")

ad_nbt <- read_csv("../useful clinical data/ad_mci_hc_neurobat_trimmed.csv")
ad_nbt <- ad_nbt %>% filter(RID %in% adni_rids)
ad_nbt$USERDATE = as.Date(ad_nbt$USERDATE, format = "%m/%d/%Y")

ad_moca <- read_csv("../useful clinical data/ad_mci_hc_moca_trimmed.csv")
ad_moca <- ad_moca %>% filter(RID %in% adni_rids)
ad_moca$USERDATE = as.Date(ad_moca$USERDATE, format = "%m/%d/%Y")

ad_faq <- read_csv("../useful clinical data/ad_mci_hc_faq_trimmed.csv")
ad_faq <- ad_faq %>% filter(RID %in% adni_rids)
ad_faq$USERDATE = as.Date(ad_faq$USERDATE, format = "%m/%d/%Y")

## nifd data -------------------------------------------------------------

nifd <- all_data[nchar(as.character(all_data$patient_original_id)) == 8, ]
nifd_ids <- unique(nifd$patient_original_id)

ftd_data <- read_csv("../useful clinical data/ftd_data.csv")
ftd_data <- ftd_data %>% filter(LONI_ID %in% nifd_ids)
ftd_data <- ftd_data[order(ftd_data$LONI_ID, decreasing = FALSE, na.last = FALSE),]
ftd_data$CLINICAL_LINKDATE <- as.Date(ftd_data$CLINICAL_LINKDATE, format = "%m/%d/%Y")

# Extracting adni data ----------------------------------------------------
## Extracting adni mmse scores --------------------------------------------------

mmse_scores = vector()

for (i in 1:length(adni_rids)) {
  
  # MMSE data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_mmse %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and mmse studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest MMSE date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the mmse score as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[1, ]
      }
      
      # Pick the closest mmse score and assign it to the imaging date
      mmse_score = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "MMSCORE"]
      
      mmse_scores = append(mmse_scores, mmse_score)
      
    } else {
      
      mmse_scores = append(mmse_scores, NA)
      
    }
  }
}

mmse_scores = unlist(mmse_scores)

## Extracting adni cdr scores --------------------------------------------------

cdr_scores = vector()

for (i in 1:length(adni_rids)) {
  
  # Clinical data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_cdr %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and clinical studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest clinical date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the clinical score as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      # Pick the closest clinical score and assign it to the imaging date
      cdr_score = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "CDGLOBAL"]
      
      cdr_scores = append(cdr_scores, cdr_score)
      
    } else {
      
      cdr_scores = append(cdr_scores, NA)
      
    }
  }
}

cdr_scores = unlist(cdr_scores)

## Extracting adni demographic information --------------------------------------------------

genders = vector()
dobs = vector()
educs = vector()

for (i in 1:length(adni_rids)) {
  
  # Clinical data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_demog %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and clinical studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest clinical date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the clinical score as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[1, ]
      }
      
      # Pick the closest clinical data and assign it to the imaging date
      gender = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "PTGENDER"]
      dob = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "DOB"]
      educ = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "PTEDUCAT"]
      
      genders = append(genders, gender)
      dobs = append(dobs, dob)
      educs = append(educs, educ)
      
    } else {
      
      genders = append(genders, NA)
      dobs = append(dobs, NA)
      educs = append(educs, NA)
      
    }
  }
}

genders = unlist(genders)
dobs = unlist(dobs)
educs = unlist(educs)

#write.csv(demog_data, "./ad_mci_hc_demog_2.csv")

## Extracting adni gds scores --------------------------------------------------

gds_scores = vector()

for (i in 1:length(adni_rids)) {
  
  # Clinical data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_gds %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and clinical studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest clinical date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the clinical data point as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      #(RID = 6912 - 79th, two tests 21 days apart)
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[1, ]
      }
      
      
      # Pick the closest GDS score and assign it to the imaging date
      gds_score = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "GDTOTAL"]
      
      gds_scores = append(gds_scores, gds_score)
      
    } else {
      
      gds_scores = append(gds_scores, NA)
      
    }
  }
}

gds_scores <- unlist(gds_scores)

## Extracting adni neurobat information --------------------------------------------------

dfws = vector()
dbws = vector()
bnts = vector()

for (i in 1:length(adni_rids)) {
  
  # Clinical data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_nbt %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and mmse studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest clinical date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the clinical score as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[1, ]
      }
      
      # Pick the closest clinical data and assign it to the imaging date
      tryCatch(dfw  <-  temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "DSPANFOR"], 
               warning = function(w) {
                 print(w)
                 print(paste(i, "&", j))
               })
      dbw = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "DSPANBAC"]
      bnt = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "BNTTOTAL"]
      
      dfws = append(dfws, dfw)
      dbws = append(dbws, dbw)
      bnts = append(bnts, bnt)
      
    } else {
      
      dfws = append(dfws, NA)
      dbws = append(dbws, NA)
      bnts = append(bnts, NA)
      
    }
  }
}

dfws = unlist(dfws)
dbws = unlist(dbws)
bnts = unlist(bnts)

## Extracting adni moca scores --------------------------------------------------

moca_scores = vector()

for (i in 1:length(adni_rids)) {
  
  # Clinical data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_moca %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and clinical studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest clinical date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the clinical data point as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[1, ]
      }
      
      # Pick the closest clinical data and assign it to the imaging date
      moca_score = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FFLUENCY"]
      
      moca_scores = append(moca_scores, moca_score)
      
    } else {
      
      moca_scores = append(moca_scores, NA)
      
    }
  }
}

moca_scores <- unlist(moca_scores)

## Extracting adni faq information --------------------------------------------------

bills = vector()
taxes = vector()
shopping = vector()
games = vector()
stove = vector()
meal = vector()
events = vector()
payattn = vector()
remdate = vector()
travel = vector()
total = vector()

for (i in 1:length(adni_rids)) {
  
  # Clinical data belonging to the patient of interest
  id = adni_rids[i]
  temp_data = ad_faq %>% filter(RID == id) %>% filter(!is.na(USERDATE))
  temp_date = temp_data$USERDATE
  
  # Imaging data belonging to the patient of interest
  index_data = adni %>% filter(RID == id) %>% filter(!is.na(date_of_study))
  index_date = index_data$date_of_study
  
  # Computing pairwise time differences between the imaging and mmse studies
  time_table = crossing(index_date, temp_date)
  time_table$diff = abs(difftime(time_table$index_date, time_table$temp_date, units = "days"))
  
  # Find the closest clinical date for each imaging date
  for (j in 1:length(index_date)){
    
    date_of_interest = index_date[j]
    date_of_interest_subset = time_table %>% filter(index_date == date_of_interest)
    time_diffs = date_of_interest_subset$diff
    
    # If there are no close dates, keep the clinical score as NA
    if (any(time_diffs < 365)) {
      
      associated_date_pair = date_of_interest_subset[time_diffs < 365, ]
      # If there are several close dates, pick the closest ones
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[associated_date_pair$diff == min(associated_date_pair$diff), ]
      }
      
      if (nrow(associated_date_pair) > 1) {
        associated_date_pair <- associated_date_pair[1, ]
      }
      
      # Pick the closest clinical data and assign it to the imaging date
      tryCatch(bls  <-  temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQFINAN"], 
               warning = function(w) {
                 print(w)
                 print(paste(i, "&", j))
               })
      txs = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQFORM"]
      shp = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQSHOP"]
      gms = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQGAME"]
      bvg = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQBEVG"]
      mel = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQMEAL"]
      evt = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQEVENT"]
      atn = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQTV"]
      rem = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQREM"]
      tvl = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQTRAVL"]
      ttl = temp_data[temp_data$USERDATE == associated_date_pair$temp_date, "FAQTOTAL"]
      
      bills = append(bills, bls)
      taxes = append(taxes, txs)
      shopping = append(shopping, shp)
      games = append(games, gms)
      stove = append(stove, bvg)
      meal = append(meal, mel)
      events = append(events, evt)
      payattn = append(payattn, atn)
      remdate = append(remdate, rem)
      travel = append(travel, tvl)
      total = append(total, ttl)
      
    } else {
      
      bills = append(bills, NA)
      taxes = append(taxes, NA)
      shopping = append(shopping, NA)
      games = append(games, NA)
      stove = append(stove, NA)
      meal = append(meal, NA)
      events = append(events, NA)
      payattn = append(payattn, NA)
      remdate = append(remdate, NA)
      travel = append(travel, NA)
      total = append(total, NA)
      
    }
  }
}

bills = unlist(bills)
taxes = unlist(taxes)
shopping = unlist(shopping)
games = unlist(games)
stove = unlist(stove)
meal = unlist(meal)
events = unlist(events)
payattn = unlist(payattn)
remdate = unlist(remdate)
travel = unlist(travel)
total = unlist(total)

## Final adni dataset -----------------------------------------------------------

all_adni_extracted_data <- cbind(class = adni$class, 
                            patient_original_id = adni$patient_original_id, 
                            RID = adni$RID,
                            demog_dos = as.character(adni$date_of_study),
                            date_of_birth = dobs,
                            biological_sex = genders, 
                            education = educs,
                            mmse_total_score = mmse_scores,
                            cdr_total_score = cdr_scores,
                            digit_span_fw = dfws, 
                            digit_span_bw = dbws, 
                            boston_naming_test = bnts,
                            letter_verbal_fluency_test = moca_scores,
                            gds_score = gds_scores,
                            faq_bills = bills,
                            faq_taxes = taxes,
                            faq_shopping = shopping,
                            faq_games = games,
                            faq_stove = stove,
                            faq_mealprep = meal,
                            faq_events = events,
                            faq_payattn = payattn,
                            faq_remdates = remdate,
                            faq_travel = travel,
                            faq_total = total)

write.csv(all_adni_extracted_data, "./adni_all_extracted_data_2.csv")


# Extracting nifd data ----------------------------------------------------

final_nifd_data = data.frame()

for (i in 1:length(nifd_ids)) {
  id = nifd_ids[i]
  nifd_subset = nifd %>% filter(patient_original_id == id)
  nifd_subset_dates = c(nifd_subset$date_of_study, NA)
  ftd_data_subset = ftd_data %>% filter(LONI_ID == id) %>% filter(CLINICAL_LINKDATE %in%  nifd_subset_dates)
  final_nifd_data = rbind(final_nifd_data, ftd_data_subset)
}
dim(final_nifd_data)

all_nifd_extracted_data = cbind(class = final_nifd_data$DX, 
                                patient_original_id = final_nifd_data$LONI_ID, 
                                demog_dos = as.character(final_nifd_data$CLINICAL_LINKDATE),
                                date_of_birth = final_nifd_data$DOB,
                                biological_sex = final_nifd_data$GENDER, 
                                education = final_nifd_data$EDUCATION,
                                mmse_total_score = final_nifd_data$MMSE_TOT,
                                cdr_total_score = final_nifd_data$CDR_TOT,
                                digit_span_fw = final_nifd_data$DIGITFW, 
                                digit_span_bw = final_nifd_data$DIGITBW, 
                                boston_naming_test = final_nifd_data$BNTCORR,
                                letter_verbal_fluency_test = final_nifd_data$DCORR,
                                gds_score = final_nifd_data$GDS15TO,
                                faq_bills = final_nifd_data$FAQ_BILLS,
                                faq_taxes = final_nifd_data$FAQ_TAXES,
                                faq_shopping = final_nifd_data$FAQ_SHOPPING,
                                faq_games = final_nifd_data$FAQ_GAMES,
                                faq_stove = final_nifd_data$FAQ_STOVE,
                                faq_mealprep = final_nifd_data$FAQ_MEALPREP,
                                faq_events = final_nifd_data$FAQ_EVENTS,
                                faq_payattn = final_nifd_data$FAQ_PAYATTN,
                                faq_remdates = final_nifd_data$FAQ_REMDATES,
                                faq_travel = final_nifd_data$FAQ_TRAVEL,
                                faq_total = final_nifd_data$FAQ_TOTAL)


write.csv(all_nifd_extracted_data, "./nifd_all_extracted_data_2.csv")


nifd_ids[!(nifd_ids %in% unique(ftd_data$LONI_ID))]
