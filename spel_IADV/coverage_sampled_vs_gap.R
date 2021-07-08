get_ent_coverage <- function(entity_id, age_start, age_end){
  
  # 1. Query the SISAL database ------------------------------------------####
  query <- paste("SELECT site.site_id, site.site_name, site.latitude, site.longitude, site.elevation, 
entity.entity_name, entity.entity_status, entity.corresponding_current, 
sample.*, gap, hiatus, d18O_measurement, d18O_precision, 
interp_age 
FROM site JOIN entity USING(site_id) 
JOIN sample USING(entity_id) 
LEFT JOIN hiatus USING(sample_id) 
LEFT JOIN gap USING(sample_id) 
LEFT JOIN original_chronology USING(sample_id) 
LEFT JOIN d18O USING(sample_id) 
WHERE entity_id =", entity_id, "AND (interp_age BETWEEN", age_start, "AND", age_end, ");", sep = " ")
  
  dt <- dbGetQuery(mydb, query)
  
  dt <- dt[order(dt$depth_sample, dt$sample_id),]
  row.names(dt) <- NULL # reset the index
  
  #x <- which(dt$interp_age <= max_age) # select samples <= 12000
  #y <- which(is.na(dt[1:x[length(x)],"interp_age"])) # select hiatuses/gaps that are <= 12000
  #xy <- sort(c(x,y)) # order both
  #dt <- dt[xy,] # filter to just <= 12000
  
  # 2. Perform filter on table queried from (1.)--------------------------####
  # Initiate an empty table to store the data from sites selected
  
  dt_site <- data.frame(site_id = integer(), sample_coverage = double(), hiatus_gap_coverage = double())
  dt_site <- c()

  all_gap_hiatus <- all((dt[['gap']] == 'G') | (dt[['hiatus']] == 'H'))
  if (is.na(all_gap_hiatus)){
    all_gap_hiatus <- F
  }
  if (all_gap_hiatus){
    next
  }
  
  # entity has gaps and/or hiatuses
  if (any((dt[['gap']] == 'G') | (dt[['hiatus']] == 'H'), na.rm = T)){
    # group the samples 
    grp_ctr <- 1
    for (k in 1:dim(dt)[1]){
      if (is.na(dt[k,'gap']) & is.na(dt[k,'hiatus'])){
        dt$grp[k] = toString(grp_ctr)
        } else {
          dt$grp[k] = NA
          grp_ctr = grp_ctr + 1
        }
      }
    dt_sum <- subset(dt, !is.na(grp)) %>% group_by(grp) %>% summarise(coverage = diff(range(interp_age, na.rm = T)), min_age = min(interp_age, na.rm = T), max_age = max(interp_age, na.rm = T))
    coverage <- sum(dt_sum[['coverage']][is.finite(dt_sum[['coverage']])], na.rm = T)
        
    gap_length_vec <- numeric()
    for (i in 1:nrow(dt[which(is.na(dt$grp) == T),])){
      gap_min <- dt[which(is.na(dt$grp) == T)[i] - 1, "interp_age"]
      gap_max <- dt[which(is.na(dt$grp) == T)[i] + 1, "interp_age"]
      if (length(gap_min) == 0 | length(gap_max) == 0){ next }
      if (gap_max > gap_min){ gap_length <- gap_max - gap_min } else {
        gap_length <- gap_min - gap_max }
          
      gap_length_vec[i] <- gap_length
    }
    gap_length_tot <- sum(gap_length_vec)
    dt_site <- rbind(dt_site, c(entity_id, coverage,gap_length_tot))
    
    } else {
      # There are no gaps or hiatuses
      coverage <- diff(range(dt$interp_age, na.rm = T))
      dt_site <- rbind(dt_site, c(entity_id, coverage, 0))
    }
  return(dt_site)
}

      
    