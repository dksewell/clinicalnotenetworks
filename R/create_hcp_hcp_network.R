#' Compute HCP-HCP network
#' 
#' Take bipartite network information and convert to HCP-HCP network,
#' respecting certain time constraints for the edges.
#' 
#' @param bp_network tbl_graph object derived from create_bp_network_data()
#' @param time_interval vector of length 2 dttm
#' 
#' @return A tidygraph object with edge weights reflecting the number 
#' of clinical notes written by one HCP and read by another HCP.
#' 
#' @export


create_hcp_hcp_network = function(bp_network,
                                  time_interval = c(mdy("1-1-1900"),
                                                    mdy("1-1-2100"))){
  
  #--------------------------------------
  # Filter nodes by time interval
  #--------------------------------------
  bp_network %<>% 
    activate(nodes) %>% 
    filter(earliest_appearance <= time_interval[2]) 
  
  
  #--------------------------------------
  # Filter edges by time interval
  #--------------------------------------
  ## Keep all modifies, and all views within correct time window
  bp_network %<>%
    activate(edges) %>% 
    filter( (EVENT_ACTION == "Modify") | 
              ( (EVENT_ACTION == "View") & 
                  (date_time >= time_interval[1]) & 
                  (date_time <= time_interval[2])) )
  
    
  #--------------------------------------
  # Create adjacency matrix
  #--------------------------------------
  
  ## Create index of hcps 
  which_hcp = 
    which( (bp_network %>% 
              activate(nodes) %>% 
              as_tibble() %>% 
              pull(type)) == "hcp")
    
  
  ## Create bipartite adjacency matrix
  bipartite_adjmat = 
    as.matrix(bp_network)
  
  ## Use matrix multiplication to get weighted hcp-hcp adjacency matrix
  hcp_hcp_adjmat = 
    bipartite_adjmat[which_hcp,-which_hcp] %*%
    bipartite_adjmat[-which_hcp,which_hcp]
  
  
  #--------------------------------------
  # Create tbl graph object
  #--------------------------------------
  
  ## Get node data frame
  node_df = 
    bp_network %>% 
    activate(nodes) %>% 
    filter(type == "hcp") %>% 
    select(-type) %>% 
    as_tibble()
  
  ## Get edge data frame
  ### Quick way leveraging sparse matrices, but
  #     without having to import igraph or Matrix 
  #     since they don't play too nicely with argon.
  get_j_indices <- function(p) {
    j_indices <- integer()
    for (j in seq_len(length(p) - 1)) {
      count <- p[j + 1] - p[j]
      j_indices <- c(j_indices, rep(j, count))
    }
    return(j_indices)
  }
  edge_df = NULL
  try({
    edge_df = 
      tibble(from = hcp_hcp_adjmat@i + 1, # Index starts at 0
             to = get_j_indices(hcp_hcp_adjmat@p),
             weight = hcp_hcp_adjmat@x)
  },silent = TRUE)
  if(is.null(edge_df)) edge_df = tibble(from = integer(0),
                                        to = integer(0))
  
  
  ## Put them together
  hcp_hcp_network = 
    tbl_graph(nodes = node_df,
              edges = edge_df)
  
  
  return(hcp_hcp_network)
}



















