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

create_bp_network_data = function(pat_id,
                                  db,
                                  site = c("UCDavis",
                                           "UCSD",
                                           "UCLA")[1],
                                  include_med_students = FALSE,
                                  only_earliest_read = TRUE,
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
           log_event_number = 1:n()) %>%  # This is to help with several operations later
    mutate(original_note_id = NEW_DATA_OBFUS_ID) %>%  # This might not be necessary, but I'd like to store the original note ID for potential future use, possibly for the next grant on viz.
    mutate(NEW_DATA_OBFUS_ID = paste(NEW_DATA_OBFUS_ID,"01",sep="-")) # Adding a -01 for the multiple authorship code





  #--------------------------------------
  # Handle notes with view as first entry
  #--------------------------------------
  # Note: This needs to go before the med student portion.  To illustrate why,
  #       consider the case where a note has its first entry as a view with the note
  #       author as a med student.  Then we need there to be a modify by that med
  #       student in the data in order to later transform the first non-med student
  #       view into a modify (which happens in the med student part of the code).
  first_entry_is_view =
    alogs %>%
    group_by(NEW_DATA_OBFUS_ID) %>%
    # arrange(date_time) %>% # Not needed because we already arranged by this.
    filter(row_number() == 1) %>%
    ungroup() %>%
    filter(EVENT_ACTION == "View")

  if( nrow(first_entry_is_view) > 0 ){
    # For each note with first entry as a view, create a new row to be
    #   added to alogs that includes a modify of the note author.  Nothing
    #   else is known, so any other choices are arbitrary.  Here's my
    #   arbitrary choice: Set the date to be 1 week prior to earliest access
    #   log date.  Everything else is set to be NA.
    for(i in 1:nrow(first_entry_is_view)){
      new_row_final =
        new_row_temp =
        alogs %>%
        filter(NEW_DATA_OBFUS_ID == first_entry_is_view$NEW_DATA_OBFUS_ID[i]) %>%
        arrange(date_time) %>%
        filter(row_number() == 1) %>%
        mutate(date_time = min(alogs$date_time) - 7 * 24 * 60^2,
               ACCESS_USER_OBFUS_ID = NOTE_AUTHOR_OBFUS_ID,
               ACCESS_USER_PROV_TYPE = NOTE_AUTHOR_PROV_TYPE,
               ACCESS_USER_CLINICIAN_TITLE = NOTE_AUTHOR_CLINICIAN_TITLE,
               ACCESS_USER_PROV_SPECIALTY = NOTE_AUTHOR_PROV_SPECIALTY,
               ACCESS_USER_PROVIDER_GENDER = NOTE_AUTHOR_PROVIDER_GENDER)
      new_row_final[] = NA
      new_row_final %<>%
        mutate(date_time = min(alogs$date_time) - 7 * 24 * 60^2,
               ACCESS_USER_OBFUS_ID = new_row_temp$NOTE_AUTHOR_OBFUS_ID,
               ACCESS_USER_PROV_TYPE = new_row_temp$NOTE_AUTHOR_PROV_TYPE,
               ACCESS_USER_CLINICIAN_TITLE = new_row_temp$NOTE_AUTHOR_CLINICIAN_TITLE,
               ACCESS_USER_PROV_SPECIALTY = new_row_temp$NOTE_AUTHOR_PROV_SPECIALTY,
               ACCESS_USER_PROVIDER_GENDER = new_row_temp$NOTE_AUTHOR_PROVIDER_GENDER,
               log_event_number = max(alogs$log_event_number) + 1)

      alogs %<>%
        bind_rows(new_row_final)
    }

    rm(new_row_final,new_row_temp)
  }




  #--------------------------------------
  # Deal with med students
  #--------------------------------------
  if(!include_med_students){

    to_be_removed = NULL

    # We'll eventually just delete all med students.  This is fine
    #   for views, but not necessarily modifies.
    med_student_authored_notes =
      alogs %>%
      filter(ACCESS_USER_PROV_TYPE == ".STUDENT: MEDICAL",
             EVENT_ACTION == "Modify")
    ## If such notes exist, change authorship.
    if(nrow(med_student_authored_notes) > 0){
      to_be_added_to_alogs = alogs[0,]

      for(j in 1:nrow(med_student_authored_notes)){
        ### First, check to see if there are any views.
        note_views =
          alogs %>%
          filter(NEW_DATA_OBFUS_ID == med_student_authored_notes$NEW_DATA_OBFUS_ID[j],
                 EVENT_ACTION == "View")

        ### Second, check for an additional modify.
        second_modify =
          alogs %>%
          filter(ACCESS_USER_PROV_TYPE != ".STUDENT: MEDICAL",
                 EVENT_ACTION == "Modify",
                 NEW_DATA_OBFUS_ID == med_student_authored_notes$NEW_DATA_OBFUS_ID[j]) %>%  # Yikes, forgot this in the previous version!
          select(ACCESS_USER_OBFUS_ID,
                 ACCESS_USER_PROV_TYPE,
                 ACCESS_USER_CLINICIAN_TITLE,
                 ACCESS_USER_PROV_SPECIALTY,
                 ACCESS_USER_PROVIDER_GENDER)

        #### If there are no other modifies nor views, this is all moot.
        #     This is because we effectively have a note with no author
        #     and no viewers.
        if( (nrow(note_views) == 0) & (nrow(second_modify) == 0) ){
          to_be_removed =
            c(to_be_removed,
              med_student_authored_notes$NEW_DATA_OBFUS_ID[j])
        }else{

          ### If there is a second modify by a non-medical student
          #     (and all other students should not be included in the db)
          #     then set that person as the new author.
          if(nrow(second_modify) > 0){
            new_author =
              second_modify[1,]
          }else{

            ### If there's not an additional modify, authorship goes to
            #     first physician view, else first view if no physician view.
            #     Add an additional row to show a link
            #     FROM healthcare provider TO note
            note_views %<>%
              mutate(phys_author = grepl("PHYSICIAN",ACCESS_USER_PROV_TYPE))
            if(any(note_views$phys_author)) note_views %<>% filter(phys_author)
            new_author = note_views[1,]
            new_author$EVENT_ACTION = "Modify"
            to_be_added_to_alogs %<>%
              bind_rows(new_author %>%
                          select(-phys_author))
          }
        }
      }

      # Add new rows onto alogs.
      to_be_added_to_alogs %<>%
        distinct() %>%
        mutate(log_event_number =
                 max(alogs$log_event_number) + 1:n())
      alogs %<>%
        bind_rows(to_be_added_to_alogs)

    }


    # OK, now delete all med student actions.
    alogs %<>%
      filter(ACCESS_USER_PROV_TYPE != ".STUDENT: MEDICAL")
    # And delete all notes with no non-student views nor modifies
    if(!is.null(to_be_removed)) alogs %<>% filter(!(NEW_DATA_OBFUS_ID %in% unique(to_be_removed)))


  }else{
    stop("If med students are to be included, then code needs to be written here to do so correctly.")
  }



  #--------------------------------------
  # Handle multiple authorship
  #--------------------------------------
  # My logic can be represented in the following example.  Suppose we have:
  # |---Note ID suffix---|---ACCESSUSEROBFUSID---|---EVENTACTION---|
  # |---------1----------|-----------A-----------|---M-------------|
  # |---------1----------|-----------B-----------|---V-------------|
  # |---------1----------|-----------B-----------|---M-------------|
  # |---------1----------|-----------B-----------|---V-------------|
  # |---------1----------|-----------B-----------|---M-------------|
  # |---------1----------|-----------D-----------|---V-------------|
  # |---------1----------|-----------E-----------|---M-------------|
  # Then the second modify will add the following rows to alogs:
  # |---------1----------|-----------B-----------|---V-------------|
  # |---------2----------|-----------B-----------|---M-------------|
  # |---------2----------|-----------C-----------|---V-------------|
  # |---------2----------|-----------D-----------|---V-------------|
  # |---------2----------|-----------E-----------|---V-------------|
  # and will remove the follow row from alogs
  # |---------1----------|-----------B-----------|---M-------------|
  # The third modify, since it is not a new author, only removes the row from alogs
  # |---------1----------|-----------B-----------|---M-------------|
  # The fourth modify will add the following rows to alogs:
  # |---------1----------|-----------E-----------|---V-------------|
  # and will remove the following row from alogs
  # |---------1----------|-----------E-----------|---M-------------|
  
  
  if(verbose) cat("\n--- Handling clinical notes with multiple authorship ---\n")
  ## Find complex notes with multiple note authors
  multiple_note_authors =
    alogs %>%
    group_by(NEW_DATA_OBFUS_ID) %>%
    filter(EVENT_ACTION == "Modify") %>%
    select(ACCESS_USER_OBFUS_ID) %>%
    distinct() %>%
    summarize(n_unique_authors = n()) %>%
    filter(n_unique_authors > 1)

  if(nrow(multiple_note_authors) > 0){

    ## Create objects to keep track of changes.  alogs will be modified at the end.
    to_be_added_to_alogs = alogs[0,]
    to_be_removed = NULL
    
    ## Create helper function to modify the DATA_OBFUS_ID
    augment_id = function(string,n){
      paste(substr(string,1,nchar(string) - 2),
            sprintf("%02d", as.numeric(substr(string,nchar(string) - 1,nchar(string))) + n),
            sep = "")
    }
    

    ## Loop through complex notes to handle multiple authorship
    for(note_id in multiple_note_authors$NEW_DATA_OBFUS_ID){

      ### Get access logs just for this note
      smaller_alogs =
        alogs %>%
        filter(NEW_DATA_OBFUS_ID == note_id)

      ### Get index of modifies
      modify_index = which(smaller_alogs$EVENT_ACTION == "Modify")

      ### Create vector of current authors (just the first one to start with.  It will expand as we go.)
      current_authors = smaller_alogs$ACCESS_USER_OBFUS_ID[modify_index[1]]

      ### Loop through these events and either add or edit the alogs
      for(i in 2:length(modify_index)){ # Start at 2 because we already used arrange(date_time)

        if( smaller_alogs$ACCESS_USER_OBFUS_ID[i] %in% current_authors ){
          #### If authorship doesn't change, remove this observation.
          to_be_removed =
            c(to_be_removed,
              smaller_alogs$log_event_number[modify_index[i]])
          # (I know this is inefficient coding, because it's in both if{} and else{}, 
          #   but conceptually it's easier to read like this)
        }else{
          #### If new author, do the following steps:
          ##### Change modify to view for original note id.
          #       Note that this effectively adds and subtracts a row.  I.e., 
          #       it removes this modify for the original note id, and adds 
          #       a new view for the original note id.
          alogs$EVENT_ACTION[which(alogs$log_event_number == 
                                     smaller_alogs$log_event_number[modify_index[i]])] = 
            "View"
          ##### Add modify for new note (increase suffix on name)
          to_be_added_to_alogs %<>%
            bind_rows(
              smaller_alogs[modify_index[i],] %>% 
                mutate(NEW_DATA_OBFUS_ID = 
                         augment_id(note_id,i - 1),
                       log_event_number = 
                         max(alogs$log_event_number) + 
                         nrow(to_be_added_to_alogs) + 1)
            )
          ##### Add views for new note.  (Only if there are subsequent 
          #       events after the last modify, of course.)
          #       Also, don't include events from the same author.
          if(modify_index[i] < nrow(smaller_alogs)){
            to_be_added_to_alogs %<>%
              bind_rows(
                smaller_alogs %>% 
                  filter(row_number() > modify_index[i],
                         ACCESS_USER_OBFUS_ID != smaller_alogs$ACCESS_USER_OBFUS_ID[modify_index[i]]) %>% # If this HCP is creating the new note, they don't need to be recorded as a viewer.
                  mutate(EVENT_ACTION = "View",
                         NEW_DATA_OBFUS_ID = 
                           augment_id(note_id,i - 1),
                         log_event_number = 
                           max(alogs$log_event_number) + 
                           nrow(to_be_added_to_alogs) + 1:n())
              )
          }
          
          #### Finally, add this author to the authorship vector.
          current_authors = 
            c(current_authors,
              smaller_alogs$ACCESS_USER_OBFUS_ID[modify_index[i]])
        }
        
      }# End loop through events in complex note
      
    }# End loop through complex notes

    ## Add additional Views to alogs
    alogs %<>%
      bind_rows(distinct(to_be_added_to_alogs))
    ## Remove certain lines from alogs
    alogs %<>%
      filter(!(log_event_number %in% to_be_removed))

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
  node_df_hcps =
    alogs %>% 
    group_by(ACCESS_USER_OBFUS_ID) %>% 
    summarize(earliest_appearance = 
                min(date_time)) %>% 
    rename(name = ACCESS_USER_OBFUS_ID) %>% 
    mutate(type = "hcp")

  ## Now combine the notes and hcps
  node_df =
    bind_rows(node_df_notes,
              node_df_hcps)


  #--------------------------------------
  # Create edge df
  #--------------------------------------

  ## Get edges which are modifies, one author column at a time.
  ### Arrange alogs by datetime since we have added rows.
  alogs %<>%
    arrange(date_time)
  ### Get unique authors
  edge_df =
    #### Get modify events
    alogs %>%
    filter(EVENT_ACTION == "Modify") %>%
    group_by(ACCESS_USER_OBFUS_ID,
             NEW_DATA_OBFUS_ID) %>% 
    filter(row_number() == 1) %>%  # I don't think this is necessary as these secondary modifies by the same author ought to have been removed above, but just in case...
    ungroup() %>% 
    select(ACCESS_USER_OBFUS_ID,
           NEW_DATA_OBFUS_ID,
           EVENT_ACTION,
           original_note_id,
           date_time) %>%
    rename(from = ACCESS_USER_OBFUS_ID,
           to = NEW_DATA_OBFUS_ID,
           !!ifelse(only_earliest_read,"earliest_date_time","date_time") := date_time)
  
  ## Get edges which are views
  edge_df %<>%
    bind_rows(
      alogs %>%
        filter(EVENT_ACTION == "View") %>%
        group_by(ACCESS_USER_OBFUS_ID,
                 NEW_DATA_OBFUS_ID) %>% 
        filter(row_number() <= ifelse(only_earliest_read,1,n())) %>% # Just get earliest event if desired.
        ungroup() %>% 
        select(NEW_DATA_OBFUS_ID,
               ACCESS_USER_OBFUS_ID,
               EVENT_ACTION,
               original_note_id,
               date_time) %>%
        rename(from = NEW_DATA_OBFUS_ID,
               to = ACCESS_USER_OBFUS_ID,
               !!ifelse(only_earliest_read,"earliest_date_time","date_time") := date_time)
    )
  

  #--------------------------------------
  # Combine to get tbl graph
  #--------------------------------------

  patient_network_data =
    tbl_graph(nodes = node_df,
              edges = edge_df)


  return(patient_network_data)
}# End main function
