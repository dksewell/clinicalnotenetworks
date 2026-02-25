#' Create network data for a patient
#'
#' This function creates a bipartite network for HCPs and
#' clinical notes.
#'
#' @details
#' NOTES:
#' \itemize{
#'  \item The curated access logs data are sufficiently small for
#'  UCDavis at least that this could be done faster by
#'  collecting all the access logs first.  However, for
#'  memory considerations I'm still doing the filtering on
#'  the patient ID on the SQLite database and then collecting.
#'  \item If the first entry in the data is a View, set the author to be
#'  the first Modify.  If no modify exists, use note author id.  The new 
#'  modify row is set to be one year prior to the patient's first 
#'  CANCER_DIAGNOSIS_DATE
#'  \item For Non-UCSD sites, change 
#'  \itemize{
#'    \item med student prov type for sites other than UCSD!
#'    \item values to look for in this_notes_alogs$ACCESS_USER_PROV_TYPE 
#'    and ACCESS_USER_CLINICIAN_TITLE for med student
#'  }
#'  \item Multiple author situations are handled as follows:
#'  \itemize{
#'    \item A new author leads to a change in the note ID going forward.
#'    \item The new note ID has multiple authors.  Previous authors
#'    carry forward.
#'    \item A new View is created with the secondary author as the
#'    reader, and the original authors as the authors.
#'  }
#' }
#'
#' @param pat_id character giving the patient id
#' @param db SQLite connection.
#' @param site Which site should we pull data for?  Should be
#' a character, namely one of "UCDavis", "UCSD", or "UCLA".
#' @param include_med_students if FALSE, med students' events are not kept.
#' @param only_earliest_read If TRUE, duplicated Views (same note, same viewer)
#' are removed and only the earliest is kept.
#' @param verbose logical.
#'
#' @return A `tidygraph` object.
#'
#' The **node tibble** contains the following:
#' \itemize{
#'   \item \code{name}: Character corresponding to the access logs and patient data.
#'   \item \code{earliest_appearance}: dttm
#'   \item \code{type}: Character, either `"note"` or `"hcp"`.
#' }
#'
#' The **edge tibble** contains the following:
#' \itemize{
#'   \item \code{from}: Integer referring back to the node tibble.
#'   \item \code{to}: Integer referring back to the node tibble.
#'   \item \code{EVENT_ACTION}: Character, either `"Modify"` or `"View"`.
#'   \item \code{date_time}: dttm
#' }
#'
#'
#' @import lubridate
#' @import dplyr
#' @import magrittr
#' @import tidygraph
#' @import RSQLite
#' @import DBI
#'
#' @export

create_bp_network_data = function(data_pat,
                                  db,
                                  site = c("UCDavis",
                                           "UCSD",
                                           "UCLA")[2],
                                  include_med_students = FALSE,
                                  verbose = TRUE){
  
  if(nrow(data_pat) > 1) stop("data_pat must have only one row")
  
  #--------------------------------------
  # Get censoring dates for the given patient
  #--------------------------------------
  
  data_pat = 
    data_pat |> 
    mutate(
      LEFT_CENSOR_DATE = as_date(LEFT_CENSOR_DATE),
      RIGHT_CENSOR_DATE = as_date(RIGHT_CENSOR_DATE),
      STUDY_ENTRY_DATE = as_date(STUDY_ENTRY_DATE),
      CANCER_DIAGNOSIS_DATE = as_date(CANCER_DIAGNOSIS_DATE))
  
  
  #--------------------------------------
  # Collect access logs for given patient
  #--------------------------------------
  if(verbose) cat("\n--- Collecting patient's access logs ---\n")
  pat_id = 
    data_pat$PAT_OBFUS_ID[1]
  
  alogs =
    db |>
    tbl(paste0("access_logs_",site)) |>
    filter(PAT_OBFUS_ID == pat_id) |>  # PAT_OBFUS_ID needs to be of same format (chr vs. numeric) as patient level data
    collect() |>
    mutate(date_time = as_datetime(ACCESS_TIME)) |> # Convert access_time_char to date
    arrange(date_time) |>
    mutate(ACCESS_USER_OBFUS_ID = as.character(ACCESS_USER_OBFUS_ID), # This will help merge with note ids later.
           log_event_number = row_number()) |>  # This is to help with several operations later
    select(date_time,
           log_event_number,
           PAT_OBFUS_ID,
           NEW_DATA_OBFUS_ID,
           EVENT_ACTION,
           ACCESS_USER_OBFUS_ID,
           ACCESS_USER_PROV_TYPE,
           ACCESS_USER_CLINICIAN_TITLE,
           ACCESS_USER_PROV_SPECIALTY,
           ACCESS_USER_PROVIDER_GENDER,
           NOTE_AUTHOR_OBFUS_ID,
           NOTE_AUTHOR_PROV_TYPE,
           NOTE_AUTHOR_CLINICIAN_TITLE,
           NOTE_AUTHOR_PROV_SPECIALTY,
           NOTE_AUTHOR_PROVIDER_GENDER)
  
  
  
  
  #--------------------------------------
  # Handle notes with view as first entry
  #--------------------------------------
  # Note: This needs to go before the med student portion.  To illustrate why,
  #       consider the case where a note has its first entry as a view with the note
  #       author as a med student.  Then we need there to be a modify by that med
  #       student in the data in order to later transform the first non-med student
  #       view into a modify (which happens in the med student part of the code).
  first_entry_is_view =
    alogs |>
    group_by(NEW_DATA_OBFUS_ID) |>
    # arrange(date_time) |> # Not needed because we already arranged by this.
    filter(row_number() == 1) |>
    ungroup() |>
    filter(EVENT_ACTION == "View")

  if( nrow(first_entry_is_view) > 0 ){
    if(verbose) cat("\n--- Handling clinical notes with first entry being a View ---\n")
    
    for(i in 1:nrow(first_entry_is_view)){
      this_notes_alogs = 
        alogs |> 
        filter(NEW_DATA_OBFUS_ID == first_entry_is_view$NEW_DATA_OBFUS_ID[i])
      
      if(any(this_notes_alogs$EVENT_ACTION == "Modify")){
        # Remove views prior to modify.  We simply don't know who wrote it.
        rows_to_remove =
          this_notes_alogs |>
          filter(cumsum(EVENT_ACTION == "Modify") == 0) |> 
          pull(log_event_number)
        alogs = 
          alogs |> 
          filter(!(log_event_number %in% rows_to_remove))
        
        # Changed this.  I disagree with this logic.
        # # If there is a modify event, set the earliest modify
        # #   to be the author.
        # 
        # earliest_modify = 
        #   this_notes_alogs |> 
        #   filter(EVENT_ACTION == "Modify") |> 
        #   filter(row_number() == 1)
        # alogs = 
        #   alogs |>
        #   add_row(
        #     date_time = left_censor_date,
        #     log_event_number = nrow(alogs) + i,
        #     PAT_OBFUS_ID = earliest_modify$PAT_OBFUS_ID,
        #     NEW_DATA_OBFUS_ID = earliest_modify$NEW_DATA_OBFUS_ID,
        #     EVENT_ACTION = "Modify",
        #     ACCESS_USER_OBFUS_ID = earliest_modify$ACCESS_USER_OBFUS_ID,
        #     ACCESS_USER_PROV_TYPE = 
        #       earliest_modify$ACCESS_USER_PROV_TYPE,
        #     ACCESS_USER_CLINICIAN_TITLE = 
        #       earliest_modify$ACCESS_USER_CLINICIAN_TITLE,
        #     ACCESS_USER_PROV_SPECIALTY = 
        #       earliest_modify$ACCESS_USER_PROV_SPECIALTY,
        #     ACCESS_USER_PROVIDER_GENDER = 
        #       earliest_modify$ACCESS_USER_PROVIDER_GENDER,
        #     NOTE_AUTHOR_OBFUS_ID = 
        #       earliest_modify$ACCESS_USER_OBFUS_ID,
        #     NOTE_AUTHOR_PROV_TYPE = 
        #       earliest_modify$ACCESS_USER_PROV_TYPE,
        #     NOTE_AUTHOR_CLINICIAN_TITLE = 
        #       earliest_modify$ACCESS_USER_CLINICIAN_TITLE,
        #     NOTE_AUTHOR_PROV_SPECIALTY = 
        #       earliest_modify$ACCESS_USER_PROV_SPECIALTY,
        #     NOTE_AUTHOR_PROVIDER_GENDER = 
        #       earliest_modify$ACCESS_USER_PROVIDER_GENDER
        #   )
      }else{
        # If there is not a modify, set the note author as the author
        
        alogs = 
          alogs |>
          add_row(
            date_time = data_pat$LEFT_CENSOR_DATE,
            log_event_number = nrow(alogs) + i,
            PAT_OBFUS_ID = first_entry_is_view$PAT_OBFUS_ID[i],
            ACCESS_USER_OBFUS_ID = first_entry_is_view$ACCESS_USER_OBFUS_ID[i],
            EVENT_ACTION = "Modify",
            NEW_DATA_OBFUS_ID = 
              first_entry_is_view$NEW_DATA_OBFUS_ID[i],
            ACCESS_USER_PROV_TYPE = 
              first_entry_is_view$NOTE_AUTHOR_PROV_TYPE[i],
            ACCESS_USER_CLINICIAN_TITLE = 
              first_entry_is_view$NOTE_AUTHOR_CLINICIAN_TITLE[i],
            ACCESS_USER_PROV_SPECIALTY = 
              first_entry_is_view$NOTE_AUTHOR_PROV_SPECIALTY[i],
            ACCESS_USER_PROVIDER_GENDER = 
              first_entry_is_view$NOTE_AUTHOR_PROVIDER_GENDER[i],
            NOTE_AUTHOR_PROV_TYPE = 
              first_entry_is_view$NOTE_AUTHOR_PROV_TYPE[i],
            NOTE_AUTHOR_CLINICIAN_TITLE = 
              first_entry_is_view$NOTE_AUTHOR_CLINICIAN_TITLE[i],
            NOTE_AUTHOR_PROV_SPECIALTY = 
              first_entry_is_view$NOTE_AUTHOR_PROV_SPECIALTY[i],
            NOTE_AUTHOR_PROVIDER_GENDER = 
              first_entry_is_view$NOTE_AUTHOR_PROVIDER_GENDER[i]
          )
      }
      
    }#End: cycling through notes having the first entry as a view
    
    alogs = 
      alogs |> 
      arrange(date_time)
    
  }#End: Handle notes having the first entry as a view




  #--------------------------------------
  # Deal with med students
  #--------------------------------------
  if(!include_med_students){
    if(verbose) cat("\n--- Handling clinical notes with medical student authorship ---\n")
    
    notes_to_be_removed = NULL

    # Find all modify events from a med student 
    med_student_authored_notes =
      alogs |>
      filter(ACCESS_USER_PROV_TYPE == "Medical Student", # Probably will need to change for sites other than UCSD
             EVENT_ACTION == "Modify")
    
    ## If such entries exist, change authorship.
    if(nrow(med_student_authored_notes) > 0){
      
      for(j in 1:nrow(med_student_authored_notes)){
        # Look at this note's entries after the med student modified it, nad
        #   remove other med student entries.
        this_notes_alogs = 
          alogs |> 
          filter(NEW_DATA_OBFUS_ID == med_student_authored_notes$NEW_DATA_OBFUS_ID[j]) |> 
          filter(date_time >= med_student_authored_notes$date_time[j],
                 ACCESS_USER_PROV_TYPE != "Medical Student")
        
        # If there is a physician involved, only consider those.  Else, look at other prov types.
        if(any(this_notes_alogs$ACCESS_USER_PROV_TYPE == "Attending Physician")){
          this_notes_alogs = 
            this_notes_alogs |> 
            filter(ACCESS_USER_PROV_TYPE == "Attending Physician")
        }else{
          if(any(this_notes_alogs$ACCESS_USER_CLINICIAN_TITLE %in% 
                 c("MD", "MD, PhD"))){
            this_notes_alogs = 
              this_notes_alogs |> 
              filter(ACCESS_USER_CLINICIAN_TITLE %in% 
                       c("MD", "MD, PhD"))
          }
        }
        
        # If there is a subsequent modify, use that.  Else, use next View
        if(any(this_notes_alogs$EVENT_ACTION == "Modify")){
          this_notes_alogs = 
            this_notes_alogs |> 
            filter(EVENT_ACTION == "Modify")
        }
        this_notes_alogs = 
          this_notes_alogs |> 
          filter(row_number() == 1)
        
        # Overwrite authorship information
        alogs_row_index = 
          which(alogs$log_event_number ==
                  med_student_authored_notes$log_event_number[j])
        
        alogs$ACCESS_USER_OBFUS_ID[alogs_row_index] = 
          this_notes_alogs$ACCESS_USER_OBFUS_ID
        alogs$ACCESS_USER_PROV_TYPE[alogs_row_index] =
          this_notes_alogs$ACCESS_USER_PROV_TYPE
        alogs$ACCESS_USER_CLINICIAN_TITLE[alogs_row_index] = 
          this_notes_alogs$ACCESS_USER_CLINICIAN_TITLE
        alogs$ACCESS_USER_PROV_SPECIALTY[alogs_row_index] = 
          this_notes_alogs$ACCESS_USER_PROV_SPECIALTY
        alogs$ACCESS_USER_PROVIDER_GENDER[alogs_row_index] = 
          this_notes_alogs$ACCESS_USER_PROVIDER_GENDER
        # alogs$NOTE_AUTHOR_OBFUS_ID # I don't think we actually care about these.
        # alogs$NOTE_AUTHOR_PROV_TYPE
        # alogs$NOTE_AUTHOR_CLINICIAN_TITLE
        # alogs$NOTE_AUTHOR_PROV_SPECIALTY
        # alogs$NOTE_AUTHOR_PROVIDER_GENDER
        
        # Remove second modify.
        alogs = 
          alogs[-which(alogs$log_event_number == this_notes_alogs$log_event_number[1]),]
      }
      
    }#End: Handling med student authorships
    
    # Remove any med student views (authorships will already be 
    #   overwritten by this point)
    alogs = 
      alogs |> 
      filter(ACCESS_USER_PROV_TYPE != "Medical Student")
    
  }#End: Handling med student entries
  
  
  #--------------------------------------
  # Create node df
  #--------------------------------------
  if(verbose) cat("\n--- Creating tidygraph object ---\n")
  
  # Create node_df for the hcps
  ## Get patient's HCP involvement times
  node_df_hcps = 
    db |> 
    tbl(paste0("involvement_times_",tolower(site))) |> 
    filter(PAT_OBFUS_ID == pat_id) |> 
    collect() |> 
    rename(name = ACCESS_USER_OBFUS_ID) |> 
    select(-PAT_OBFUS_ID) |> 
    mutate(type = "hcp",
           date_start = as_date(date_start),
           date_end = as_date(date_end))
  
  
  ## Create node_df for the notes
  node_df_notes =
    alogs |>
    group_by(NEW_DATA_OBFUS_ID) |>
    summarize(earliest_appearance =
                min(date_time)) |>
    rename(name = NEW_DATA_OBFUS_ID) |>
    mutate(type = "note")
    

  ## Now combine the notes and hcps
  node_df =
    bind_rows(node_df_notes,
              node_df_hcps)


  #--------------------------------------
  # Create edge df
  #--------------------------------------

  ## Get edges which are modifies, one author column at a time.
  alogs =
    alogs |> 
    arrange(date_time)
  ### Get unique authors
  edge_df =
    #### Get modify events
    alogs |>
    filter(EVENT_ACTION == "Modify") |>
    group_by(ACCESS_USER_OBFUS_ID,
             NEW_DATA_OBFUS_ID) |> 
    filter(row_number() == 1) |>  # This might require some discussion.
    ungroup() |> 
    select(ACCESS_USER_OBFUS_ID,
           NEW_DATA_OBFUS_ID,
           EVENT_ACTION,
           date_time) |>
    rename(from = ACCESS_USER_OBFUS_ID,
           to = NEW_DATA_OBFUS_ID)
  
  ## Get edges which are views
  edge_df =
    edge_df |> 
    bind_rows(
      alogs |>
        filter(EVENT_ACTION == "View") |>
        select(NEW_DATA_OBFUS_ID,
               ACCESS_USER_OBFUS_ID,
               EVENT_ACTION,
               date_time) |>
        rename(from = NEW_DATA_OBFUS_ID,
               to = ACCESS_USER_OBFUS_ID)
    )
  

  #--------------------------------------
  # Combine to get tbl graph
  #--------------------------------------

  patient_network_data =
    tbl_graph(nodes = node_df,
              edges = edge_df)


  return(patient_network_data)
}# End main function
