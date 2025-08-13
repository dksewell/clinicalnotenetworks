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
#'  \item I have filtered out telephone encounters
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

create_bp_network_data = function(pat_id,
                                  db,
                                  site = c("UCDavis",
                                           "UCSD",
                                           "UCLA")[1],
                                  verbose = TRUE){
  
  #--------------------------------------
  # Collect access logs for given patient
  #--------------------------------------
  if(verbose) cat("\n--- Collecting patient's access logs ---\n")
  alogs = 
    db %>% 
    tbl(paste0("access_logs_",site)) %>% 
    filter(PAT_OBFUS_ID == pat_id, # PAT_OBFUS_ID needs to be of same format (chr vs. numeric) as patient level data
           NOTE_TYPE != "Telephone Encounter") %>%  
    collect() %>% 
    mutate(date_time = dmy_hms(ACCESS_TIME_CHAR)) %>% # Convert access_time_char to date
    arrange(date_time) %>% 
    mutate(ACCESS_USER_OBFUS_ID = as.character(ACCESS_USER_OBFUS_ID), # This will help merge with note ids later.
           log_event_number = 1:n()) # This is to help with several operations later
  
  # We will be adding rows (views), so we need to keep track of the
  #   number of rows of alogs.
  max_log_event_number = nrow(alogs)
  to_be_added_to_alogs = alogs[0,] # This creates an empty tibble of the correct form
  
  
  #--------------------------------------
  # Handle multiple authorship
  #--------------------------------------
  if(verbose) cat("\n--- Handling clinical notes with multiple authorship ---\n")
  ## Find complex notes
  complex_notes_ids =
    alogs %>% 
    group_by(NEW_DATA_OBFUS_ID) %>% 
    filter(EVENT_ACTION == "Modify") %>% 
    select(ACCESS_USER_OBFUS_ID) %>% 
    distinct() %>% 
    summarize(n_unique_authors = n()) %>% 
    filter(n_unique_authors > 1)
  
  ## Create new chr columns of alogs for multiple authors
  if(nrow(complex_notes_ids) > 0){
    for(j in 2:max(complex_notes_ids$n_unique_authors)){
      alogs[[paste0("NOTE_AUTHOR_OBFUS_ID_",j)]] = character(nrow(alogs))
    }
  
    ## Loop through complex notes to handle multiple authorship
    for(note_id in complex_notes_ids$NEW_DATA_OBFUS_ID){
      
      ### Get access logs just for this note
      smaller_alogs =
        alogs %>% 
        filter(NEW_DATA_OBFUS_ID == note_id)
      
      ### Loop through these events and either add or edit the alogs
      for(i in 2:nrow(smaller_alogs)){ # Start at 2 because we already used arrange(date_time)
        
        #### Only act for modifies
        if( smaller_alogs$EVENT_ACTION[i] != "Modify" ) next
        
        #### Only act if the authorship changes
        if( smaller_alogs$ACCESS_USER_OBFUS_ID[i] %in% 
            (smaller_alogs[i,] %>% 
             select(contains("NOTE_AUTHOR_OBFUS_ID")) %>% 
             unlist() %>% 
             unique()) ) next
        # If we haven't run next yet, we have a new modifier!
        
        #### Add a row for a view on the original note ID.
        ####  Change the event action to view, and change the log_event_number
        to_be_added_to_alogs %<>%
          bind_rows(smaller_alogs[i,] %>% 
                      mutate(EVENT_ACTION = "View",
                             log_event_number = 
                               max_log_event_number + 1))
        #### Update the maximum log_event_number for next time
        max_log_event_number = max_log_event_number + 1
        
        #### Create new note id moving forward
        smaller_alogs$NEW_DATA_OBFUS_ID[i:nrow(smaller_alogs)] = 
          paste0(smaller_alogs$NEW_DATA_OBFUS_ID[1],"_",i) # This won't be sequential (1,2,3,...) but it will be unique
        #### Add new author to authorship list moving forward
        ##### Get the earliest column which is currently blank
        author_columns = 
          smaller_alogs[i,] %>% 
          select(contains("NOTE_AUTHOR_OBFUS_ID"))
        earliest_blank_author_col = min(which(unlist(author_columns) == ""))
        ##### Add author to this column for all subsequent entries
        smaller_alogs[[colnames(author_columns)[earliest_blank_author_col]]][i:nrow(smaller_alogs)] =
          smaller_alogs$ACCESS_USER_OBFUS_ID[i]
        
      }# End loop through events in complex note
      
      ### Edit alogs based on smaller_alogs
      alogs[smaller_alogs$log_event_number,] = 
        smaller_alogs
      
      
    }# End loop through complex notes
    
    ## Add additional Views to alogs
    alogs %<>%
      bind_rows(to_be_added_to_alogs)
    
  }
  
  #--------------------------------------
  # Create node df
  #--------------------------------------
  if(verbose) cat("\n--- Creating tidygraph object ---\n")
  
  ## Create node_df for the notes
  node_df_notes = 
    alogs %>% 
    group_by(NEW_DATA_OBFUS_ID) %>% 
    summarize(earliest_appearance = 
                min(date_time)) %>% 
    rename(name = NEW_DATA_OBFUS_ID) %>% 
    mutate(type = "note")
  
  ## Create node_df for the hcps
  ### Get unique IDs
  node_df_hcps =
    tibble(name = 
             alogs %>% 
             select(ACCESS_USER_OBFUS_ID,
                    contains("NOTE_AUTHOR_OBFUS_ID")) %>% 
             unlist() %>%
             unique() %>% 
             setdiff(""),
           earliest_appearance = mdy_hms("1-1-2100_00:00:00"),
           type = "hcp")
  ### Get earliest appearance in the access logs
  for(i in 1:nrow(node_df_hcps)){
    node_df_hcps$earliest_appearance[i] = 
      alogs %>%
      rowwise() %>% 
      filter(any(c_across(c("ACCESS_USER_OBFUS_ID",starts_with("NOTE_AUTHOR_OBFUS_ID"))) == node_df_hcps$name[i])) %>% 
      ungroup() %>% 
      pull(date_time) %>% 
      min()
  }
  
  ## Now combine the notes and hcps
  node_df = 
    bind_rows(node_df_notes,
              node_df_hcps)
  
  
  #--------------------------------------
  # Create edge df
  #--------------------------------------
  
  ## Get edges which are modifies, one author column at a time.
  ### Get unique authors
  edge_df =
    #### Get modify events
    alogs %>% 
    filter(EVENT_ACTION == "Modify") %>% 
    select(ACCESS_USER_OBFUS_ID,
           NEW_DATA_OBFUS_ID,
           EVENT_ACTION,
           date_time) %>% 
    rename(from = ACCESS_USER_OBFUS_ID,
           to = NEW_DATA_OBFUS_ID) %>% 
    #### Some events may not have the authorship entry if,
    #       e.g., the note was authored prior to the data
    #       collection time period
    bind_rows(
      alogs %>% 
        select(NOTE_AUTHOR_OBFUS_ID,
               NEW_DATA_OBFUS_ID,
               EVENT_ACTION,
               date_time) %>% 
        rename(from = NOTE_AUTHOR_OBFUS_ID,
               to = NEW_DATA_OBFUS_ID)
    )
  #### Continue with last task, but account for complex notes
  if(nrow(complex_notes_ids) > 0){
    for(j in 2:max(complex_notes_ids$n_unique_authors)){
      edge_df %<>%
        bind_rows(
          alogs %>% 
            select(!!paste0("NOTE_AUTHOR_OBFUS_ID","_",j),
                   NEW_DATA_OBFUS_ID,
                   EVENT_ACTION,
                   date_time) %>% 
            filter(get(paste0("NOTE_AUTHOR_OBFUS_ID","_",j)) != "",
                   EVENT_ACTION == "View") %>% 
            rename(from = !!paste0("NOTE_AUTHOR_OBFUS_ID","_",j),
                   to = NEW_DATA_OBFUS_ID)
        )
    }
  }
  #### Keep only unique authorship edges.  First (earliest) 
  #     entry will be kept.
  edge_df %<>%
    distinct(to,from,.keep_all = TRUE)
  
  ## Get edges which are views
  edge_df %<>%
    bind_rows(
      alogs %>% 
      filter(EVENT_ACTION == "View") %>% 
      select(NEW_DATA_OBFUS_ID,
             ACCESS_USER_OBFUS_ID,
             EVENT_ACTION,
             date_time) %>% 
      rename(from = NEW_DATA_OBFUS_ID,
             to = ACCESS_USER_OBFUS_ID)
    )
  # Note: We do **NOT** want to use distinct on the views!
  
  
  #--------------------------------------
  # Combine to get tbl graph
  #--------------------------------------
  
  patient_network_data = 
    tbl_graph(nodes = node_df,
              edges = edge_df)
  
  
  return(patient_network_data)
}# End main function
