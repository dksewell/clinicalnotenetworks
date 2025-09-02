#' Feature extraction for both-sides regression modeling
#' 
#' Find patient level variables, note level features, and group readership
#' indicators.
#' 
#' @param pat_id character giving the patient id 
#' @param bp_network tbl_graph object derived from create_bp_network_data()
#' @param time_before_dx integer.  Number of days prior to date of diagnosis.
#' @param end_time dttm
#' @param inpatient_prop_thresh Numeric between 0 and 1.  What proportion of 
#' access logs should be during an inpatient stay to consider the HCP on an 
#' inpatient team?  Default is 0.75 based on looking at histograms for UCDavis.
#' @param db SQLite connection.
#' @param site Which site should we pull data for?  Should be 
#' a character, namely one of "UCDavis", "UCSD", or "UCLA".
#' 
#' @importFrom readr parse_number
#' @export


feature_extraction = function(pat_id,
                              bp_network,
                              time_before_dx = 30,
                              end_time,
                              inpatient_prop_thresh = 0.75,
                              db,
                              site = c("UCDavis",
                                       "UCSD",
                                       "UCLA")[1]){
  
  #--------------------------------------
  # Get baseline features
  #--------------------------------------
  
  
  # Note: OTHER_THERAPY_SUMMARY_DATE_STAR is too rare to include.
  pat_features = 
    db %>% 
    tbl(paste0("patient_data_",
               site)) %>% 
    filter(PAT_OBFUS_ID == pat_id) %>% 
    collect() %>% 
    mutate(dx_date = ymd("1970-01-01") + TRUE_DX_DATE) %>% 
    mutate(end_time = end_time) %>% 
    # Compute whether treatment has been initiated
    mutate(surgery = 
             factor(isTRUE(ymd("1899-12-30") + as.integer(SURGERY_FIRST_COURSE_DATE_OF_R) <= end_time),
                    levels = c(T,F)),
           radiation = 
             factor(isTRUE(ymd("1899-12-30") + as.integer(RADIATION_DATE_STARTED_REP_1) <= end_time),
                    levels = c(T,F)),
           chemo = 
             factor(isTRUE(ymd("1899-12-30") + as.integer(RADIATION_DATE_STARTED_REP_1) <= end_time),
                    levels = c(T,F)),
           hormone = 
             factor(isTRUE(ymd("1899-12-30") + as.integer(HORMONE_SUMMARY_DATE_STARTED) <= end_time),
                    levels = c(T,F)),
           immunotherapy = 
             factor(isTRUE(ymd("1899-12-30") + as.integer(IMMUNOTHERAPY_SUMMARY_DATE_STAR) <= end_time),
                    levels = c(T,F))) %>% 
    # Get cancer site and stage
    mutate(cancer_site = 
             case_match(SITE_GROUP,
                        "BREAST" ~ "breast",
                        c("COLON","RECTUM & RECTOSIGMOID") ~ "colon",
                        "LUNG/BRONCHUS-NON SM CELL" ~ "lung") %>% 
             factor(levels = c("breast","colon","lung")),
           cancer_stage = 
             readr::parse_number(TNM_MIXED_STAGE) %>% 
             factor(levels = 1:4)) %>% 
    # Get comorbidities
    select(-CANCER_LEUK,-CANCER_LYMPH,-CANCER_METS,-CANCER_NSITU,-CANCER_SOLID) %>% 
    mutate(across(AIDS:WGHTLOSS,
                  ~ factor(isTRUE(dmy(.x) <= end_time),
                           levels = c(T,F)))) %>% 
    # Select all the above (+ demographics)
    select(PAT_OBFUS_ID,
           dx_date,
           end_time,
           AGE_AT_DIAGNOSIS,
           SEX,
           RACE,
           SPANISH_HISPANIC_ORIGIN,
           surgery,
           radiation,
           chemo,
           hormone,
           immunotherapy,
           cancer_site,
           cancer_stage,
           AIDS:WGHTLOSS) %>% 
    rename_with(tolower)
  
  
  
  #--------------------------------------
  # Filter edges by time interval
  #--------------------------------------
  
  bp_network %<>%
    activate(edges) %>% 
    filter(date_time <= end_time,
           date_time >= pat_features$dx_date - time_before_dx)
  
  
  #--------------------------------------
  # Filter nodes by time interval
  #--------------------------------------
  
  bp_network %<>% 
    activate(nodes) %>% 
    filter(earliest_appearance <= end_time) %>% 
    filter(centrality_degree(mode = "all") > 0) # This is because there may be some HCPs from too far back not involved in patient's cancer care
  
  
  
  #--------------------------------------
  # Match up MTS group
  #--------------------------------------
  
  # Get unique mts teams
  mts_teams = 
    db %>% 
    tbl(paste0("mts_lookup_",site)) %>% 
    select(mts_component_group) %>% 
    distinct() %>% 
    collect()
  mts_teams %<>%
    bind_rows(tibble(mts_component_group = 
                       paste(mts_teams$mts_component_group,
                             "inpatient",
                             sep="_"))
    ) %>% 
    pull(mts_component_group)
  
  # Collect the mts lookup table for the patient
  lu = 
    db %>% 
    tbl(paste0("mts_lookup_",site)) %>% 
    filter(PAT_OBFUS_ID == pat_id) %>% 
    collect()
  
  # Convert inpatient teams
  lu %<>%
    mutate(mts_component_group = 
             ifelse(proportion_logs_inpatient < inpatient_prop_thresh,
                    mts_component_group,
                    paste(mts_component_group,"inpatient",sep = "_")))
  
  # Join them up to node df
  bp_network %<>%
    activate(nodes) %>% 
    left_join(
      lu %>% 
        select(ACCESS_USER_OBFUS_ID,
               mts_component_group),
      by = join_by(name == ACCESS_USER_OBFUS_ID)
    )
  
  # Finally, join them up to edge df
  bp_network %<>%
    activate(edges) %>% 
    mutate(from_name = .N()$name[from],
           to_name = .N()$name[to]) %>% 
    left_join(
      lu %>% 
        select(ACCESS_USER_OBFUS_ID,
               mts_component_group) %>% 
        rename(from_mts_component_group = mts_component_group),
      by = join_by(from_name == ACCESS_USER_OBFUS_ID)
    ) %>% 
    left_join(
      lu %>% 
        select(ACCESS_USER_OBFUS_ID,
               mts_component_group) %>% 
        rename(to_mts_component_group = mts_component_group),
      by = join_by(to_name == ACCESS_USER_OBFUS_ID)
    ) %>% 
    select(-from_name,-to_name)
  
  
  #--------------------------------------
  # Get note importance
  #--------------------------------------
  
  # Extract time since dx as a covariate
  note_importance = 
    bp_network %>% 
    activate(edges) %>% 
    filter(EVENT_ACTION == "Modify") %>% 
    mutate(time_from_dx = as.numeric(as.Date(date_time) - pat_features$dx_date)) %>% 
    select(time_from_dx,from_mts_component_group) %>% 
    as.data.frame()
  
  # # Match up mts group of author to a new variable, one for each mts group
  # for(j in mts_teams){
  #   note_importance %<>%
  #     mutate(!!j := from_mts_component_group == j)
  # }
  # Nevermind, I don't think I want to do it this way, because we'll want to
  #   construct a design matrix from this using the from_mts_comp._group as
  #   a factor variable.  We'll also want an intercept.
  
  # (Clean up columns via select() later)
  
  # NOTE: I've left essentially double the number of MTS groups here on purpose.  
  #   We can always condense down to the original groups, and add a single 
  #   variable for inpatient/outpatient designation based on this should we
  #   wish to later for parsimony.
  
  
  #--------------------------------------
  # Get note readings
  #--------------------------------------
  
  # This who_read_it matrix is to have \phi subtracted off.
  who_read_it = 
    matrix(0.0,
           nrow(note_importance),
           NROW(mts_teams),
           dimnames = list(NULL,
                           mts_teams))
  for(i in 1:nrow(note_importance)){
    note_indices = 
      bp_network %>% 
      activate(edges) %>% 
      filter(from == note_importance$to[i]) %>% 
      as.data.frame()
    
    if(nrow(note_indices) > 0){
      who_read_it[i,unique(note_indices$to_mts_component_group)] = 2.0
    }
  }
  
  
  
  list(PAT_OBFUS_ID = pat_features$pat_obfus_id,
       end_time = end_time,
       patient_features = pat_features,
       note_importance = 
         note_importance %>% 
         select(-from,-to) %>% 
         mutate(from_mts_component_group = 
                  factor(from_mts_component_group,
                         levels = mts_teams)),
       who_read_it = who_read_it) %>% 
    return()
  
  
}