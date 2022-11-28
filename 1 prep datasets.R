library(dplyr)
library(data.table)
library(progress)
library(data.table)
library(stringr)
library(tidyr)

paths_base = dir('results', full.names = TRUE)

output = lapply(paths_base, function(name_dataset_this) {
  print(name_dataset_this)
  
  filenames_all_dataset_this = dir(path=name_dataset_this, recursive = TRUE, full.names = TRUE)
  
  # make a table of all the runs
  df_fn = rbindlist(lapply(filenames_all_dataset_this, function(fn) {
    splits = strsplit(gsub("capped_2", "capped2", basename(fn)),split="_")[[1]]
    
    if (length(splits)<5) # if we have the adjacenc/assemblage case
    {
      splits_raw = splits
      splits = rep(NA, 7)
      splits[1] = splits_raw[1]
      splits[2] = splits_raw[2]
      splits[7] = splits_raw[3]
    }
    
    result = data.frame(
      dataset=basename(name_dataset_this),
      filename = fn,
      epsilon = gsub("1.0e.5","0.00001",gsub("-",".",gsub("eps","",splits[1]))) %>% as.numeric,
      capped = splits[2]=="capped2",
      cost_add = gsub("-",".",gsub("add","",splits[3])) %>% as.numeric,
      cost_delete = gsub("-",".",gsub("del","",splits[4])) %>% as.numeric,
      cost_environment = gsub("-",".",gsub("temp","",splits[5])) %>% as.numeric,
      cost_wait = gsub("\\.txt","",gsub("-",".",gsub("wait","",splits[6]))) %>% as.numeric,
      
      type =  gsub("^a$","astar",gsub("\\.csv","",gsub("\\.txt","",splits[7])))
      )
    })) %>%
    group_by(dataset, epsilon, cost_add, cost_delete, cost_environment, cost_wait, capped)
  
  df_fn_aa = df_fn %>% 
    filter(type %in% c("adjacency","assemblage")) %>%
    ungroup
  
  df_fn = df_fn %>% 
    filter(!(type %in% c("adjacency","assemblage"))) %>%
    mutate(filename_astar = filename[type=="astar"], 
           filename_summary = filename[type=="summary"], 
           filename_transitions = filename[type=="transitions"]) %>%
    slice(1) %>%
    select(-filename,-type)
  
  # get stats for each row
  df_stats = rbindlist(lapply(df_fn$filename_summary, function(fn) {
    df_summary = read.csv(fn,sep="\n", header=FALSE)
    stats = sapply(strsplit(df_summary[c(13,15,19:37),],split=": "), function(x) { z = x[2]; names(z) = x[1]; return(z)  }) %>% t %>% as.data.frame
    return(stats)
    })) %>% mutate_at(3:ncol(.), as.numeric)
  
  df_all = cbind(df_fn, df_stats) %>% as.data.frame

  return(list(df_all=df_all,df_fn_aa=df_fn_aa))
})


df_all = rbindlist(lapply(output, function(x) {x$df_all}))
df_aa = rbindlist(lapply(output, function(x) {x$df_fn_aa}))


# get transition statistics
df_transitions = rbindlist(lapply(1:nrow(df_all), function(i)
{
  df_transitions_this = read.csv(df_all$filename_transitions[i])
  intermediate_visits_q95 = quantile(df_transitions_this$intermediate_visits,0.95)
  in_degree_q95 = quantile(df_transitions_this$in_degree,0.95)
  out_degree_q95 = quantile(df_transitions_this$out_degree,0.95)
  
  return(data.frame(intermediate_visits_q95, in_degree_q95, out_degree_q95))
}))

df_correlation_r2 = rbindlist(lapply(1:nrow(df_all), function(i)
{
  df_transitions_this = read.csv(df_all$filename_transitions[i])
  r2.intermediate_visits_in_degree = summary(lm(in_degree~intermediate_visits,data=df_transitions_this))$r.squared
  
  print(i/nrow(df_all))
  return(data.frame(r2.intermediate_visits_in_degree))
}))

paths_counts = lapply(1:nrow(df_all), function(i) {
  print(i/nrow(df_all))
  df_astar_this = fread(df_all$filename_astar[i]) # fread is faster than read.csv
  
  MAX_COUNT = 20 # arbitrary value
  
  vals_this = df_astar_this$a_star_ops
  vals_this[is.na(vals_this)] = MAX_COUNT # pick a big number
  
  table_this = tabulate(vals_this)
  result = data.frame(t(table_this))
  names(result)[MAX_COUNT]="NOT_REACHABLE"
  return(list(counts=result, df=df_astar_this))
})

df_astar_fractions = rbindlist(lapply(1:length(paths_counts), function(i) { paths_counts[[i]]$counts / sum(paths_counts[[i]]$counts) }))
names(df_astar_fractions) = paste("astar_length",names(df_astar_fractions),sep="_")

df_astar_paths_nontrivial = lapply(1:length(paths_counts), function(i) { paths_counts[[i]]$df })
rm(paths_counts)



# get abundances
df_aa_ss = df_aa %>% filter(type=='assemblage')
df_abundances = rbindlist(lapply(1:nrow(df_aa_ss), function(i)
{
    df_assemblages_this = read.csv(df_aa_ss$filename[i])
    
    abundance_grand_mean = mean(df_assemblages_this$abundance_mean,na.rm=T)
    abundance_grand_sd = sd(df_assemblages_this$abundance_mean,na.rm=T)
    
    return(data.frame(abundance_grand_mean, abundance_grand_sd))
})) %>% 
  cbind(df_aa_ss %>% select(dataset)) %>%
  unique






df_all_joined = cbind(df_all, df_transitions, df_correlation_r2, df_astar_fractions)
df_all_joined_with_abundance = df_all_joined %>% left_join(df_abundances,by='dataset')


if (!file.exists('outputs'))
{
  dir.create('outputs')
}

write.csv(df_all_joined_with_abundance, file='outputs/df_all_joined.csv', row.names=FALSE)
write.csv(df_aa, file='outputs/df_aa.csv', row.names=FALSE)



# join in the degree info for the paths
df_astar_paths_nontrivial_with_degree = lapply(1:nrow(df_all_joined), function(i) {
  print(i/nrow(df_all_joined))
  
  transitions_this = read.csv(df_all_joined$filename_transitions[i])
  paths_this = df_astar_paths_nontrivial[[i]]
  paths_this_final = paths_this %>% left_join(transitions_this %>% 
                             rename(start=node) %>% 
                             select(start, initial_out_degree=out_degree)) %>% 
                                    left_join(transitions_this %>% 
                  rename(goal=node) %>% 
                  select(goal, desired_in_degree=in_degree))
  
  return(paths_this_final)
})




# save all the paths too
df_astar_paths_nontrivial_condensed = lapply(df_astar_paths_nontrivial_with_degree, function(df) {
  df %>% select(start:full_path, a_star_cost, a_star_ops, initial_out_degree, desired_in_degree)  
})


# this takes about 60 seconds
saveRDS(df_astar_paths_nontrivial_condensed, file='outputs/df_astar_paths_nontrivial_condensed.Rdata',compress = TRUE)
