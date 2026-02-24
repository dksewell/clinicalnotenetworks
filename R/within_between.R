#' Compute within/between communication measures
#' 
#' @details
#' For each team at each time point, compute
#' 1. Team size: Number of team members
#' 2. Information production: Number of notes written by team members
#' 3. Within-team Information Reach: Number of note reads (excluding re-reading) 
#' from notes written by other team members
#' 4. Within-team Information share: Proportion of reads that are actualized
#' from notes written by someone else on the same team
#' 5a. Missed opportunities: Number of potential note reads - actual note reads
#' 5b. Missed opportunities: Number of completely unread notes
#' 
#' and for each other team
#' 6. Between-team Information Reach: Number of note reads (excluding re-reading) 
#' from notes written by someone on the source team
#' 7. Between-team Information Share: Proportion of reads that are actualized
#' from notes written by someone else on the source team
#' 8a. Missed opportunities: Number of potential note reads - actual note reads
#' 8b. Missed opportunities: Number of completely unread notes
#' 
#' 
#' @import DBI
#' @import RSQLite
#' @import lubridate
#' @import tidygraph
#' @import tibble 
#' @import dplyr
#' @importFrom utils data
#' @export



within_between = function(id_pat,
                          db,
                          bp_path = "/Shared/ehr_networks/R01-SMART_cancer_care/patient_networks/ucsd/bipartite/",
                          site = "UCSD",
                          date_start,
                          date_end){
  
  # Pull in data
  
  ## Read in involvement data
  data_involvement = 
    db |> 
    tbl(paste0("involvement_times_",site)) |> 
    filter(PAT_OBFUS_ID == id_pat) |> 
    collect() |> 
    mutate(date_start = as_date(date_start),
           date_end = as_date(date_end)) |> 
    mutate(date_interval = interval(date_start,date_end))
  
  
  ## Read in mts lookup table
  data_mts = 
    db |> 
    tbl(paste0("mts_lookup_",site)) |> 
    filter(PAT_OBFUS_ID == id_pat) |> 
    collect() |> 
    filter(mts != "Exclude")
  
  
  ## Match MTS group to data_involvement
  #   (Note that this will also remove anyone without an MTS assignment.)
  data_involvement = 
    data_involvement |> 
    left_join(
      data_mts |> 
        select(ACCESS_USER_OBFUS_ID,
               mts), # No need for PAT_OBFUS_ID, as it has already been filtered to this patient
      by = "ACCESS_USER_OBFUS_ID"
    ) |> 
    drop_na() |> 
    filter(mts != "Exclude") |>  # Just a check in case someone changes code above
    mutate(mts = paste(mts,setting,sep="_")) # Consider an mts group to be the combo of group and setting
  
  
  ## Read in bipartite data
  data_bipartite = 
    readRDS(paste0(bp_path,
                   "bp_data-",
                   id_pat,
                   ".RDS")) |> 
    ### Filter out nodes we don't intend to keep
    activate(nodes) |>
    filter( (type == "note") |
              (name %in% data_involvement$ACCESS_USER_OBFUS_ID) )
    
  ### Create edge df with mts groups assigned
  data_edges = 
    data_bipartite |> 
    activate(edges) |> 
    # Add ACCESS_USER_ID's
    mutate(id_from = .N()$name[from],
           id_to = .N()$name[to]) |> 
    as_tibble() |> 
    select(-from,-to) |> 
    # Add MTS group_setting
    left_join(
      data_involvement |> 
        select(ACCESS_USER_OBFUS_ID,
               mts) |> 
        rename(id_from = ACCESS_USER_OBFUS_ID,
               mts_from = mts),
      by = "id_from"
    ) |> 
    left_join(
      data_involvement |> 
        select(ACCESS_USER_OBFUS_ID,
               mts) |> 
        rename(id_to = ACCESS_USER_OBFUS_ID,
               mts_to = mts),
      by = "id_to"
    )
  
  ## Read in patient data
  data_pat = 
    db |> 
    tbl(paste0("patient_data_",site)) |> 
    filter(PAT_OBFUS_ID == id_pat) |> 
    collect()
  
  
  # Get date sequences
  if(missing(date_start)){
    date_start = 
      data_pat$STUDY_ENTRY_DATE |> 
      as_date()
  }
  if(missing(date_end)){
    date_end = 
      data_pat$RIGHT_CENSOR_DATE |> 
      as_date()
  }
  ## Lookback 4 weeks if possible
  max_lookback =
    floor(
      as.integer(date_start - 
                   as_date(data_pat$LEFT_CENSOR_DATE)) / 7
    )
  ## Create sequence of dates
  seq_week = 
    (-min(max_lookback,4)):floor(as.integer(date_end - date_start) / 7)
  seq_startdates =
    date_start + weeks(seq_week)
    
  
  # Create large tibble, one row for each day, one column for each network measure
  
  ## Initialize tibble
  network_measures = 
    tibble(start = 
             seq_startdates[-length(seq_startdates)],
           end = seq_startdates[-1] - 1,
           week = seq_week[-length(seq_week)])
  
  
  ## get unique MTS groups
  utils::data(mts,package = "clinicalnotenetworks")
  unique_mts = 
    unique(mts$mts) |> 
    setdiff("Exclude")
  unique_mts = 
    c(paste(unique_mts,"outpatient",sep="_"),
      paste(unique_mts,"inpatient",sep="_"))
  
  ## get unique paris of MTS groups
  unique_mts_pairs = 
    expand.grid(unique_mts,unique_mts) |> 
    filter(Var1 != Var2) |> 
    mutate(mts_pairing = 
             paste(Var1,Var2,sep=":")) |> 
    pull(mts_pairing)
  
  ## Add within team measures
  ### 1. Team size: Number of team members
  for(j in paste("within_size",unique_mts,sep="_")){
    network_measures[[j]] = numeric(nrow(network_measures))
  }
  
  ### 2. Information production: Number of notes written by team members
  for(j in paste("within_production",unique_mts,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 3. Within-team Information Reach: Number of note reads (excluding re-reading) 
  #' from notes written by other team members
  for(j in paste("within_reach",unique_mts,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 4. Within-team Information share: Mean/median (over notes) proportion of 
  #       team members who read a note written by someone else on the team
  #       and for each other team
  for(j in paste("within_share",unique_mts,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 5a. Missed opportunities: Number of potential note reads - actual note reads
  for(j in paste("within_missedoppsA",unique_mts,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 5b. Missed opportunities: Number of completely unread notes
  for(j in paste("within_missedoppsB",unique_mts,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 6. Between-team Information Reach: Number of note reads (excluding re-reading) 
  #       team members who read a note written by someone else on the source team
  for(j in paste("between_reach",unique_mts_pairs,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 7. Between-team Information Share: Mean/median (over notes) proportion of 
  #       team members who read a note written by someone else on the source team
  for(j in paste("between_share",unique_mts_pairs,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 8a. Missed opportunities: Number of potential note reads - actual note reads
  for(j in paste("between_missedoppsA",unique_mts_pairs,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  ### 8b. Missed opportunities: Number of completely unread notes
  for(j in paste("between_missedoppsB",unique_mts_pairs,sep="_")){
    network_measures[[j]] = NA * numeric(nrow(network_measures))
  }
  
  
  
  
  # Cycle through each week.
  
  ## Store the teams which have appeared at least once
  applicable_teams_cum = NULL
  
  for(wk in 1:nrow(network_measures)){
    ## Create time interval
    week_interval = 
      interval(network_measures$start[wk],
               network_measures$end[wk])
    
    
    ## Find applicable teams
    data_involvement_week = 
      data_involvement |> 
      filter(
        int_overlaps(date_interval,
                     week_interval)
      )
    
    ## Move on to next week if there are no applicable teams...
    if(nrow(data_involvement_week) == 0) next
    
    ## ...else record applicable teams
    applicable_teams =
      data_involvement_week |> 
      pull(mts) |> 
      unique()
    applicable_teams_cum = 
      union(applicable_teams_cum,
            applicable_teams)
    
    for(j in applicable_teams){
      ## Get applicable team members
      team_members = 
        data_involvement_week |> 
        filter(mts == j) |> 
        pull(ACCESS_USER_OBFUS_ID) |> 
        unique()
      
      ## Add within team measures
      ### 1. Team size: Number of team members
      network_measures[[paste("within_size",j,sep="_")]][wk] = 
        length(team_members)
      
      ### 2. Information production: Number of notes written by team members
      network_measures[[paste("within_production",j,sep="_")]][wk] = 
        data_edges |> 
        filter(as_date(date_time) %in% week_interval,
               mts_from == j) |> 
        pull(id_to) |> 
        unique() |> # Don't double count multiple modifies on the same note 
        length()
      
      ### 3. Within-team Information Reach: Number of note reads (excluding re-reading) 
      #       from notes written by other team members
      #       (cumulative)
      # and 
      ### 4. Within-team Information share: Mean/median (over notes) proportion of 
      #       team members who read a note written by someone else on the team
      #       and for each other team
      #       (cumulative)
      # and
      # 5a. Missed opportunities: Number of potential note reads - actual note reads
      
      ### Pull writes from team j
      from_team_j = 
        data_edges |> 
        filter(mts_from == j,
               as_date(date_time) <= int_end(week_interval)) |> 
        arrange(date_time)
      
      if(nrow(from_team_j) > 0){
        ### only keep earliest modify for each author-note combo.  
        #     This can't be run if the above code for from_team_j has 
        #     no rows.
        from_team_j = 
          from_team_j|>
          group_by(id_from,id_to) |> 
          mutate(keep = c(T,rep(F,n() - 1))) |> 
          ungroup() |> 
          filter(keep) |> 
          select(-keep)
        
        ### Keep track of the potential and actual number of reads by active team members
        potential_reads = 
          actual_reads = 
          0L
        
        ### Cycle through each writing event done by a member of team j
        for(writing_event in 1:nrow(from_team_j)){
          
          note_id = from_team_j$id_to[writing_event]
          
          potential_reads =
            potential_reads +
            # = Number of team members, excluding the note author
            length(setdiff(team_members,
                           from_team_j$id_from[writing_event]))
          
          actual_reads =
            actual_reads +
            data_edges |> 
            #### only consider subsequent actions
            filter(as_date(date_time) <= int_end(week_interval),
                   date_time > from_team_j$date_time[writing_event]) |> 
            #### remove the note author as a viable reader
            filter( (id_from != from_team_j$id_from[writing_event]) &
                      (id_to != from_team_j$id_from[writing_event]) )|> 
            #### include both subsequent views AND subsequent modifies as reads
            filter( ( (id_from %in% team_members) &
                        (id_to == note_id) ) | # Modifier for this note
                      ( (id_to %in% team_members) &
                          (id_from == note_id) ) ) |> # Viewer for this note
            #### keep unique readers
            select(id_from,
                   id_to) |> 
            unlist() |> 
            intersect(team_members) |> # removes duplicates and removes notes
            length()
          
          
        }#End: cycle through all writing events from team j
        
        ### Compute 3. w/in team reach
        network_measures[[paste("within_reach",j,sep="_")]] = 
          actual_reads
        ### Compute 4. w/in team share
        network_measures[[paste("within_share",j,sep="_")]] = 
          actual_reads / 
          potential_reads
        ### Compute 5a. w/in team missed opportunities 
        network_measures[[paste("within_missedoppsA",j,sep="_")]] = 
          actual_reads - potential_reads
        
        
        ### Compute 5b. w/in team missed opportunities 
        network_measures[[paste("within_missedoppsB",j,sep="_")]] = 
          data_edges |> 
          #### Filter to only notes written by team j
          mutate(note_id = 
                   ifelse(is.na(mts_from),
                          id_from,id_to)) |> 
          filter(note_id %in% from_team_j$id_to) |> 
          #### Did any team member (not an author) View it? (not view-by-subsequent-modify)
          group_by(note_id) |> 
          summarize(unread = 
                      ##### There has to be a team member available to read it who is not also an author
                      (length(setdiff(team_members,id_from)) > 0) &
                      ##### Look to see if there was a view by a team member
                      ( length(intersect(team_members,
                                         id_to)) == 0) ) |> 
          pull(unread) |> 
          sum()
        
          
      }#End: check to see if any writes by team j
      
      
      ## Add between team measures
      ### Look at other teams that might be sending info
      applicable_teams_other = 
        setdiff(applicable_teams_cum,j)
      if(length(applicable_teams_other) > 0){
        
        ### Cycle through other teams
        for(k in applicable_teams_other){
          
          #### Pull writes from team k
          from_team_k = 
            data_edges |> 
            filter(mts_from == k,
                   as_date(date_time) <= int_end(week_interval)) |> 
            arrange(date_time)
          
          if(nrow(from_team_k) > 0){
            #### only keep earliest modify 
            from_team_k = 
              from_team_k |> 
              group_by(id_from,id_to) |> 
              mutate(keep = c(T,rep(F,n() - 1))) |> 
              ungroup() |> 
              filter(keep) |> 
              select(-keep)
            
            #### Keep track of the potential and actual number of reads by active team members
            potential_reads = 
              actual_reads = 
              0L
            
            #### Cycle through each writing event done by a member of team k
            for(writing_event in 1:nrow(from_team_k)){
              
              note_id = from_team_k$id_to[writing_event]
              
              potential_reads =
                potential_reads +
                length(team_members)
              
              actual_reads =
                actual_reads +
                data_edges |> 
                ##### only consider subsequent actions
                filter(as_date(date_time) <= int_end(week_interval),
                       date_time > from_team_k$date_time[writing_event]) |> 
                ##### include both subsequent views AND subsequent modifies as reads
                filter( ( (id_from %in% team_members) &
                            (id_to == note_id) ) | 
                          ( (id_to %in% team_members) &
                              (id_from == note_id) ) ) |> 
                ##### keep unique readers
                select(id_from,
                       id_to) |> 
                unlist() |> 
                intersect(team_members) |> 
                length()
              
              
            }
            
            #### Compute 6. btwn team reach
            network_measures[[paste("between_reach",j,sep="_")]] = 
              actual_reads
            #### Compute 7. btwn team share
            network_measures[[paste("between_share",j,sep="_")]] = 
              actual_reads / 
              potential_reads
            ### Compute 8a. w/in team missed opportunities 
            network_measures[[paste("between_missedoppsA",j,sep="_")]] = 
              actual_reads - potential_reads
            
            
            ### Compute 8b. w/in team missed opportunities 
            network_measures[[paste("between_missedoppsB",j,sep="_")]] = 
              data_edges |> 
              #### Filter to only notes written by team j
              mutate(note_id = 
                       ifelse(is.na(mts_from),
                              id_from,id_to)) |> 
              filter(note_id %in% from_team_k$id_to) |> 
              #### Did any team member (not an author) View it? (not view-by-subsequent-modify)
              group_by(note_id) |> 
              summarize(unread = 
                          ##### There has to be a team member available to read it who is not also an author
                          (length(setdiff(team_members,id_from)) > 0) &
                          ##### Look to see if there was a view by a team member
                          ( length(intersect(team_members,
                                             id_to)) == 0) ) |> 
              pull(unread) |> 
              sum()
            
            
          }#End: check if any notes are written by team k
          
          
        }#End: cycle through other teams as sources
        
        
      }#End: check if any other applicable teams besides j
      
      
    }#End: cycling through each applicable team for this week
    
    
  }#End: cycling through all weeks
  
  
  
  return(network_measures)
}