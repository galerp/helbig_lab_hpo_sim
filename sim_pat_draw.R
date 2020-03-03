sim_pat_draw <- function(sim_score_file,task_id,permutation)  {

 num_pats <- pat_pairs[[task_id]]  
  estimate_mode <- function(x){
    d <- density(x)
    d$x[which.max(d$y)]
  }
  
  sim_score = read.csv(sim_score_file,stringsAsFactors=F) %>% select(-X)
  sim_score<- do.call(data.frame, lapply(sim_score, function(x) {
    replace(x, is.infinite(x) | is.na(x), 0)
  })
  )
  sim_score[is.na(sim_score)] <- 0
  cor_names = names(sim_score)
  names(sim_score) = cor_names
  rownames(sim_score) = cor_names
  
  r_100k = as.data.frame(matrix(nrow = permutation, ncol = 3))
  names(r_100k) = c("median","mean", "mode")
  
  pat_vect = names(sim_score)
  
  for(n in 1: nrow(r_100k)){
    IDs = sample(pat_vect, num_pats)
    sub_sim = sim_score[(rownames(sim_score) %in% IDs), (names(sim_score) %in% IDs)]
    diag(sub_sim) = 12345
    vect_scores =  unlist(sub_sim)
    vect_scores = vect_scores[-which(vect_scores == 12345)]
    
    r_100k$median[n] = median(vect_scores)
    r_100k$mean[n] = mean(vect_scores)
    r_100k$mode[n] = estimate_mode(vect_scores)
  }
  # write.csv(r_100k,paste("100k_",num_pats,"mod_epgp-egrp_v10.csv", sep = ""),row.names = F)
  return(r_100k)
}

