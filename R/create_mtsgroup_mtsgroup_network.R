#' Create MTS group-MTS group network
#' 
#' Take the bipartite information and create a network connecting
#' the MTS groups.  Connections and absences are now persistent over time.
#' 
#' @param pat_id character giving the patient id 
#' @param bp_network tbl_graph object derived from create_bp_network_data()
#' @param end_time dttm last date of network construction
#' @param inpatient_prop_thresh Numeric between 0 and 1.  What proportion of 
#' access logs should be during an inpatient stay to consider the HCP on an 
#' inpatient team?  Default is 0.75 based on looking at histograms for UCDavis.
#' @param db SQLite connection.
#' @param site Which site should we pull data for?  Should be 
#' a character, namely one of "UCDavis", "UCSD", or "UCLA".
#' 
#' @return A tidygraph object with edge weights reflecting the number 
#' of clinical notes written by one HCP and read by another HCP.
#' 
#' @export

create_mtsgroup_mtsgroup_network = function(pat_id,
                                            bp_network,
                                            end_time,
                                            inpatient_prop_thresh = 0.75,
                                            db,
                                            site = c("UCDavis",
                                                     "UCSD",
                                                     "UCLA")[1]){
  
  #--------------------------------------
  # Filter nodes by time interval
  #--------------------------------------
  
  bp_network %<>% 
    activate(nodes) %>% 
    filter(earliest_appearance <= end_time) 
  
  
  
  #--------------------------------------
  # Filter edges by time interval
  #--------------------------------------
  
  bp_network %<>%
    activate(edges) %>% 
    filter( date_time <= end_time )
  
  
  
  #--------------------------------------
  # Match up MTS group
  #--------------------------------------
  
  # First, collect the mts lookup table
  lu = 
    db %>% 
    dbReadTable(paste0("mts_lookup_",site))
  
  # Convert inpatient teams
  lu %<>%
    mutate(mts_component_group = 
             ifelse(proportion_logs_inpatient < inpatient_prop_thresh,
                    mts_component_group,
                    paste(mts_component_group,"inpatient",sep = "_")))
  
  # Filter to this patient
  lu %<>%
    filter(PAT_OBFUS_ID == pat_id)
  
  # Second, join them up to node df
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
  # Create MTS group-MTS group network
  #--------------------------------------
  
  # Create new node df
  node_df = 
    bp_network %>% 
    activate(nodes) %>% 
    filter(type == "hcp") %>% 
    as_tibble() %>% 
    select(mts_component_group) %>% 
    distinct() %>% 
    rename(name = mts_component_group)
  
  # Create new edge df
  ## Create the frame for a dense signed network
  edge_df =
    expand.grid(from = node_df$name,
                to = node_df$name) %>% 
    # filter(from != to) %>% # Actually, self-loops should be considered here.
    mutate(positive_weight = 0L,
           negative_weight = 0L)
  
  ## Grab the original edges
  orig_edges = 
    bp_network %>% 
    activate(edges) %>% 
    as_tibble() %>% 
    arrange(EVENT_ACTION,date_time,from,to)
  
  ## Cycle through the notes written to add + or - edges
  for(i in which(orig_edges$EVENT_ACTION == "Modify")){
    source_id = 
      orig_edges$from_mts_component_group[i]
    note_id = 
      orig_edges$to[i]
    positive_target_ids = 
      orig_edges$to_mts_component_group[which(orig_edges$from ==
                                                note_id)]
    negative_target_ids = 
      setdiff(node_df$name,positive_target_ids)
    
    edge_df %<>%
      mutate(positive_weight = 
               positive_weight + 
               ifelse( (from == source_id) &
                         (to %in% positive_target_ids),
                       1L, 0L),
             negative_weight = 
               negative_weight - 
               ifelse( (from == source_id) &
                         (to %in% negative_target_ids),
                       1L, 0L))
  }
  
  ## Get rid of any non-edges (only possible if a group didn't write any notes)
  edge_df %<>%
    filter( (positive_weight != 0) | (negative_weight != 0) )
  
  
  # Combine to get tidygraph
  return(
    tbl_graph(nodes = node_df,
              edges = edge_df)
  )
  
  
}
















