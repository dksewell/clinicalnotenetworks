#' Assign HCP's to a MTS Component Group
#' 
#' Take Level 1 MTS Master List, clean it, and apply it to access log data
#' to perform a hierarchical assignment algorithm.
#'
#' 
#' @param master_list A n x 2 .xlsx file with n component groups and a list of
#' provider specialties, clinician titles, and provider types
#' @param data data frame which includes as fields:
#' \itemize{
#'  \item ACCESS_USER_PROV_SPECIALTY
#'  \item ACCESS_USER_CLINICIAN_TITLE
#'  \item ACCESS_USER_PROV_TYPE
#' } 
#' 
#' @return A list with two elemments. The first is data.frame object/lookup 
#' table consisting of the patient id, access user's id, the access user's role 
#' type, title and specialty, and their final MTS component group assignment.
#' The second is a data.frame with the same variables as the lookup table, but 
#' the combinations are all lacking assignments.
#' 
#' @export


assign_mts <- function(master_list, data){
 
 # rename columns of the master_list
 names(mts) <- c("ID", "MTS", "SPEC_TITLE_TYPE")
 
 # Separate all specialties/titles/role types
 specialty_list <- list()
 
 for(i in 1:nrow(master_list)){
   
   specialty_list[[i]] <- master_list[i,3] %>%
     str_split(pattern = "; ") %>%
     as.data.frame()
   
   specialty_list[[i]]$MTS <- master_list$MTS[i]
   
   names(specialty_list[[i]]) <- c("SPEC_TITLE_TYPE", "MTS")
   
 }
 
 # convert to a data frame for easier merging
 mts_level_1 <- do.call(bind_rows, specialty_list)
 
 # filter each mapping to specialty, type and title, and merge to logs
 access_user_specialties <- mts_level_1 %>%
   filter(SPEC_TITLE_TYPE %in% logs$ACCESS_USER_PROV_SPECIALTY) %>%
   rename(ACCESS_USER_PROV_SPECIALTY = "SPEC_TITLE_TYPE", ACCESS_USER_MTS_SPECIALTY = "MTS") %>%
   distinct()
 
 
 access_logs <- access_logs %>%
   left_join(access_user_specialties, by = "ACCESS_USER_PROV_SPECIALTY")
 
 
 access_user_titles <- mts_level_1 %>%
   filter(SPEC_TITLE_TYPE %in% logs$ACCESS_USER_CLINICIAN_TITLE) %>%
   rename(ACCESS_USER_CLINICIAN_TITLE = "SPEC_TITLE_TYPE", ACCESS_USER_MTS_TITLE = "MTS") %>%
   distinct()
 
 
 access_logs <- access_logs %>%
   left_join(access_user_titles, by = "ACCESS_USER_CLINICIAN_TITLE")
 
 
 
 access_user_types <- mts_level_1 %>%
   filter(SPEC_TITLE_TYPE %in% logs$ACCESS_USER_PROV_TYPE) %>%
   rename(ACCESS_USER_PROV_TYPE = "SPEC_TITLE_TYPE", ACCESS_USER_MTS_TYPE = "MTS") %>%
   distinct()
 
 
 access_logs <- access_logs %>%
   left_join(access_user_types, by = "ACCESS_USER_PROV_TYPE")
 
 
 # Apply the hierarchical algorithm
 
 ## Initialize the variable
 access_logs$ACCESS_USER_MTS <- "Missing"
 
 ## First assign based on the presence of the specialty 
 ## (ideally this is if we're missing the specialty to begin)
 
 for(i in 1:nrow(access_logs)){
   if(is.na(access_logs$ACCESS_USER_PROV_SPECIALTY[i]) == F){
     access_logs$ACCESS_USER_MTS[i] <- 
       access_logs$ACCESS_USER_MTS_SPECIALTY[i]
   }
 }
 
 ## Now assign based on clinician title if specialty assignment is missing
 for(i in 1:nrow(access_logs)){
   if(access_logs$ACCESS_USER_MTS[i] == "Missing"){
     access_logs$ACCESS_USER_MTS[i] <- access_logs$ACCESS_USER_MTS_TITLE[i]
   }
 }
 
 ## Finally assign based on role type
 for(i in 1:nrow(access_logs)){
   if(access_logs$ACCESS_USER_MTS[i] == "Missing"){
     access_logs$ACCESS_USER_MTS[i] <- access_logs$ACCESS_USER_MTS_TYPE[i]
   }
 }
 
 
 # Construct the lookup table for all patient-HCP-type-title-specialty
 # combos with a component group assignment
 
 logs_with_teams <- access_logs %>%
   filter(ACCESS_USER_MTS != "Missing") %>%
   select(PAT_OBFUS_ID, ACCESS_USER_OBFUS_ID, ACCESS_USER_PROV_TYPE,
          ACCESS_USER_CLINICIAN_TITLE, ACCESS_USER_PROV_SPECIALTY,
          ACCESS_USER_MTS) %>%
   distinct()
 
 
 logs_without_teams <- access_logs %>%
   filter(ACCESS_USER_MTS == "Missing") %>%
   select(PAT_OBFUS_ID, ACCESS_USER_OBFUS_ID, ACCESS_USER_PROV_TYPE,
          ACCESS_USER_CLINICIAN_TITLE, ACCESS_USER_PROV_SPECIALTY,
          ACCESS_USER_MTS) %>%
   distinct()
 
 
 mts_lookups <- list(logs_with_teams, logs_without_teams)
 
 return(mts_lookups)
 
}

