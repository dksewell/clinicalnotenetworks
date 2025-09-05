#' Clean up patient data
#' 
#' Clean up cases in patient data with multiple rows.  
#' 
#' @details
#' \itemize{
#'  \item Clean up dates
#'  \item Clean up cancer site and stage variables
#'  \item Change DATE_LAST_PATIENT_CONTACT_OR_DEA to dx of next cancer site, if applicable
#'  \item Change stages to date intervals.
#' }
#' 
#' @param patient_data Collected patient_data_site from smart_cancer_care.db
#' 
#' @export



clean_patient_data = function(patient_data){
  
  # Convert date variables from characters
  patient_data %<>% 
    mutate(DATE_LAST_PATIENT_CONTACT_OR_DEA = ymd("1899-12-30") + as.integer(DATE_LAST_PATIENT_CONTACT_OR_DEA),
           TRUE_DX_DATE = ymd("1970-01-01") + as.integer(TRUE_DX_DATE))
  
  # Clean up Site and Stage variables
  patient_data %<>%
    mutate(cancer_site = 
             case_match(SITE_GROUP,
                        "BREAST" ~ "breast",
                        c("COLON","RECTUM & RECTOSIGMOID") ~ "colon",
                        "LUNG/BRONCHUS-NON SM CELL" ~ "lung") %>% 
             factor(levels = c("breast","colon","lung")),
           cancer_stage = 
             readr::parse_number(TNM_MIXED_STAGE) %>% 
             factor(levels = 1:4))
  
  
  # Arrange by patient ID, then date of dx.
  patient_data %<>% 
    arrange(PAT_OBFUS_ID,
            TRUE_DX_DATE)
  
  # Create variables for stage date ranges
  patient_data %<>%
    mutate(stage_1_range = interval(ymd("2100-01-01"),ymd("2101-01-01")),
           stage_2_range = interval(ymd("2100-01-01"),ymd("2101-01-01")),
           stage_3_range = interval(ymd("2100-01-01"),ymd("2101-01-01")),
           stage_4_range = interval(ymd("2100-01-01"),ymd("2101-01-01")))
  
  # Handle duplicate patient rows
  
  ## Get list of patient IDs that need to be handled
  duplicated_pat_ids = 
    patient_data %>% 
    filter(duplicated(PAT_OBFUS_ID)) %>% 
    distinct() %>% 
    pull(PAT_OBFUS_ID)
  
  ## Cycle through IDs and handle
  #   The general strategy is to create a single new row, and replace
  #   the existing rows with it.
  if(length(duplicated_pat_ids) > 0){
    for(i in duplicated_pat_ids){
      ### Create temporary data frame
      temp_df = 
        patient_data %>% 
        filter(PAT_OBFUS_ID == i)
      
      ### Check for second site and set the dx date as patient's new end date
      if(any(temp_df$cancer_site[-1] != temp_df$cancer_site[1])){
        temp_df$DATE_LAST_PATIENT_CONTACT_OR_DEA[1] = 
          min( temp_df$TRUE_DX_DATE[ which(temp_df$cancer_site != temp_df$cancer_site[1]) ] )
        temp_df %<>%
          filter(cancer_site == temp_df$cancer_site[1])
      }
      if(nrow(temp_df) > 1){ 
        ### Check for multiple stages.
        #### Get rid of redundant stages (somehow this happens; e.g., 16153930587).
        temp_df %<>%
          mutate(nonredundant = 
                   c(TRUE,
                     temp_df$cancer_stage[-nrow(temp_df)] ==
                       temp_df$cancer_stage[-1])) %>% 
          filter(nonredundant)
        #### Get stage ranges
        if(nrow(temp_df) > 1){
          stage_starts = rep(ymd("2100-01-01"),4)
          stage_ends = rep(ymd("2101-01-01"),4)
          
          ##### Initialize start and stop period for first stage in data
          stage_starts[as.integer(temp_df$cancer_stage[1])] = 
            temp_df$TRUE_DX_DATE[1]
          stage_ends[as.integer(temp_df$cancer_stage[1])] = 
            temp_df$DATE_LAST_PATIENT_CONTACT_OR_DEA[1]
          
          for(r in 2:nrow(temp_df)){
            ##### Update this row's stage interval
            stage_starts[as.integer(temp_df$cancer_stage[r])] =
              min(c(stage_starts[as.integer(temp_df$cancer_stage[r])],
                    temp_df$TRUE_DX_DATE[r])) # min() in case a stage is repeated, e.g., stage 2,3, then back to 2
            stage_ends[as.integer(temp_df$cancer_stage[r])] =
              temp_df$DATE_LAST_PATIENT_CONTACT_OR_DEA[r]
            
            ##### Update prior rows' stage interval
            prior_stages = 
              setdiff(unique(as.integer(temp_df$cancer_stage[1:(r-1)])),
                      as.integer(temp_df$cancer_stage[r]))
            stage_ends[prior_stages] = 
              sapply(stage_ends[prior_stages],
                     function(x) min(x,temp_df$TRUE_DX_DATE[r] - 1))
          }
          ##### update temp_df stage intervals
          temp_df$stage_1_range[1] = 
            interval(stage_starts[1],stage_ends[1])
          temp_df$stage_2_range[1] = 
            interval(stage_starts[2],stage_ends[2])
          temp_df$stage_3_range[1] = 
            interval(stage_starts[3],stage_ends[3])
          temp_df$stage_4_range[1] = 
            interval(stage_starts[4],stage_ends[4])
            
        }
      }
      
      ## Remove old rows and add new one
      patient_data %<>%
        filter(PAT_OBFUS_ID != temp_df$PAT_OBFUS_ID) %>% 
        bind_rows(temp_df[1,])
    }
  }
  
  
  patient_data
}






















