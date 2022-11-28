library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(data.table)
library(igraph)
library(colorRamps)
library(ranger)
library(pdp)
library(tidyverse)

df_all_joined = read.csv('outputs/df_all_joined.csv')
df_astar_paths_nontrivial = readRDS('outputs/df_astar_paths_nontrivial_condensed.Rdata')


######################################
if (!file.exists('figures'))
{
  dir.create('figures')
}
######################################

df_names = data.frame(name=c('Bucci','Carrara','Maynard','Maynard15-17-19-21-23','Maynard15-19-23','Venturelli'),
                       nice_name=c('Mouse gut','Protist','Ciliate','Ciliate+environment5','Ciliate+environment3','Human gut'),
                       n=c(11,11,5,5,5,12),
                       m=c(1,1,1,5,3,1)) %>%
  left_join(df_all_joined %>% 
              select(dataset, Proportion.of.candidate.states, abundance_grand_mean, abundance_grand_sd) %>% 
              unique %>%
              rename(name=dataset)) %>%
  mutate(Proportion.of.candidate.states=format(Proportion.of.candidate.states,digits=1),
         abundance_grand_mean=format(abundance_grand_mean,digits=1),
         abundance_grand_sd=format(abundance_grand_sd,digits=1)) %>%
  rename(`Abundance (grand mean)`=abundance_grand_mean,
         `Abundance (grand s.d.)`=abundance_grand_sd,
         `Proportion of states feasible+stable`=Proportion.of.candidate.states,
         `Dataset`=nice_name,
         `Number of species (n)`=n,
         `Number of environments (m)`=m) %>%
  arrange(Dataset)
write.csv(df_names, file='figures/df_names.csv',row.names=F)



################ DRAW A STATE DIAGRAM
df_aa = read.csv('outputs/df_aa.csv')


###################################################################################################################

classify_transients_internal <- function(df_adj)
{
  t_start = df_adj$start
  t_action = df_adj$a_str
  
  parts_start = gsub("\\*", " ", sapply(strsplit(t_start, "\\|"), head, 1))
  
  parts_action = rbindlist(lapply(strsplit(t_action, ")", fixed=TRUE), function(x) { 
    x_add = x[grep(x, pattern="\\+")] 
    if (length(x_add)==0)
    {
      x_str_add = ""
    }
    else
    {
      x_str_add = gsub("(+", "", x_add, fixed=TRUE)
    }
    
    x_del = x[grep(x, pattern="\\-")] 
    if (length(x_del)==0)
    {
      x_str_del = ""
    }
    else
    {
      x_str_del = gsub("(-", "", x_del, fixed=TRUE)
    }
    
    
    return(data.frame(parts_add=x_str_add, parts_del=x_str_del))
  }))
  
  df_combined = cbind(parts_start, t_action, parts_action)
  
  return(df_combined)
}

# try to get graph structure
classify_transients <- function(fn_adj, fn_ass)
{
  df_adj = read.csv(fn_adj)
  df_ass = read.csv(fn_ass)
  
  df_combined = classify_transients_internal(df_adj)
    
  states_transient_add = sapply(1:nrow(df_combined), function(i) {
    union(strsplit(df_combined$parts_start[i],split=" +")[[1]], 
          strsplit(df_combined$parts_add[i],split=" +")[[1]]) %>%
      as.numeric %>%
      sort %>%
      paste(collapse="*")
  })
  
  df_ass_combined = df_ass %>% 
    mutate(candidate=as.logical(candidate)) %>% 
    mutate(is.transient = gsub("|1","",str, fixed=TRUE) %in% states_transient_add) # this only works for single-environment visualizations!
  
  return(df_ass_combined)
}


blue2red_a <- function(n,a)
{
  cols = blue2red(n)
  alpha = c(0:9,letters[1:6])[floor(16*a)]
  cols_out = paste(cols,alpha,alpha,sep="")
  return(cols_out)
}


plot_state_diagram <- function(dataset_this, 
                               epsilon_this, 
                               capped_this,
                               cost_add, 
                               cost_delete,
                               cost_env,
                               cost_wait,
                               cex_vertex=1,
                               PLOT_SMALL, 
                               PLOT_LABELS,
                               DO_INVASION_GRAPH=FALSE,
                               DO_UNINVASION_GRAPH=FALSE)
{
  df_aa_ss = df_aa %>% filter(epsilon==epsilon_this & dataset==dataset_this & capped==capped_this)
  
  fn_adj = df_aa_ss %>% filter(type=='adjacency') %>% pull(filename)
  fn_ass = df_aa_ss %>% filter(type=='assemblage') %>% pull(filename)
  
  transitions_this = read.csv(fn_adj)
  
  if (DO_INVASION_GRAPH==TRUE)
  {
    transitions_this = transitions_this %>%
      filter(add==1 & del==0 & temp==0)
  }
  else if (DO_UNINVASION_GRAPH==TRUE)
  {
    transitions_this_uninvasion = transitions_this %>% filter(temp==0) # no env change
    env_this_uninvasion = sapply(strsplit(transitions_this_uninvasion$start, split="|", fixed=TRUE), tail, 1)
    df_combined_transients = classify_transients_internal(transitions_this_uninvasion)
    states_transient_add = sapply(1:nrow(df_combined_transients), function(i) {
      union(strsplit(df_combined_transients$parts_start[i],split=" +")[[1]], 
            strsplit(df_combined_transients$parts_add[i],split=" +")[[1]]) %>%
        as.numeric %>%
        sort %>%
        paste(collapse="*")
    }) %>%
      paste(env_this_uninvasion, sep="|")  
    
    transitions_this = transitions_this_uninvasion %>% 
      select(final, add, del, temp, a_str) %>% 
      cbind(start=states_transient_add) %>%
      select(start, final, everything()) %>%
      filter(start!=final)
  }
  
  states_this = classify_transients(fn_adj, fn_ass)
  
  costs_this = transitions_this$add * cost_add +
                transitions_this$del * cost_delete +
                transitions_this$temp * cost_env +
                cost_wait
  print(summary(costs_this))
  alpha_this = 0.5*(1.0 - 0.5 + 0.5 * costs_this / max(costs_this))
  
  g = graph_from_data_frame(transitions_this %>% select(start, final, everything()), vertices = states_this %>% select(str, everything()))
  
  if (PLOT_SMALL==TRUE)
  {
    g = g - which(as.logical(V(g)$candidate)==FALSE)
  }
  
  coords_this = data.frame(state=V(g)$state, richness=V(g)$richness) %>% 
    mutate(row_id = 1:nrow(.)) %>%
    group_by(richness) %>%
    mutate(x=1 - order(state)/(1+max(order(state)))) %>%
    ungroup %>%
    mutate(y=richness/max(richness))
  
  add_del_index = -1*(E(g)$add - E(g)$del)/(E(g)$add + E(g)$del)
  
  plot(g, 
       layout = coords_this %>% select(x,y) %>% as.matrix,
       vertex.frame.color = NA,
       vertex.label = if(PLOT_LABELS==TRUE) {gsub("|1", "", gsub('*', ' ', V(g)$name, fixed=TRUE), fixed=TRUE)} else {NA},
       vertex.shape='circle',
       vertex.label.family='sans',
       vertex.label.color='black',
       vertex.color=ifelse(as.logical(V(g)$stable) & as.logical(V(g)$feasible) ,'green',ifelse(V(g)$is.transient & !DO_INVASION_GRAPH, 'darkorange','gray')),
       vertex.size=cex_vertex*ifelse(as.logical(V(g)$stable) & as.logical(V(g)$feasible) ,1.25,0.75),#sqrt(scale(ifelse(is.nan(V(g)$abundance_mean), 0, V(g)$abundance_mean), center=FALSE)*5),
       edge.arrow.size=0.25,
       edge.label=if(PLOT_LABELS==TRUE) { E(g)$a_str} else {NA},
       edge.label.family='sans',
       edge.label.color=blue2red_a(100,alpha_this)[cut(add_del_index, breaks=seq(-1,1,length.out=100), include.lowest=TRUE)],
       edge.color=blue2red_a(100,alpha_this)[cut(add_del_index, breaks=seq(-1,1,length.out=100), include.lowest=TRUE)],
       edge.label.cex=0.5,
       edge.width=1.0/(costs_this),
       rescale=TRUE,
       edge.curved=FALSE
  )
}

# pdf(file='figures/g_state_diagram_example_bucci.pdf',width=15,height=15)
# plot_state_diagram(dataset_this='Bucci', epsilon_this=0.1, capped_this = FALSE,
#                    cost_add = 1,
#                    cost_delete = 10,
#                    cost_env = 1,
#                    cost_wait = 0.1,
#                    PLOT_SMALL = FALSE, PLOT_LABELS = FALSE)
# dev.off()

png(file='figures/g_state_diagram_all.png',width=10*3*100,height=10*2*100,res=100)
par(mfrow=c(2,3))
par(mar=c(0,0,4,0))
lapply(1:length(df_names$name), function(i) 
  {
  print(i)
  plot_state_diagram(dataset_this=df_names$name[i], epsilon_this=0.1, capped_this = FALSE,
                     cost_add = 1,
                     cost_delete = 1,
                     cost_env = 1,
                     cost_wait = 0.1,
                     cex_vertex = 15.0 / df_names$`Number of species (n)`[i],
                     PLOT_SMALL = FALSE, PLOT_LABELS = FALSE)
  title(sprintf("(%s) %s", letters[i], df_names$Dataset[i]),adj=0,cex.main=3,line=0)
})
dev.off()

png(file='figures/g_state_diagram_all_epsilon.png',width=10*6*100,height=10*3*100,res=100)
par(mfrow=c(3,6))
par(mar=c(0,0,4,0))
epsilons=c(1e-5,1e-3,1e-1)
lapply(epsilons, function(e)
{
  print(e)
  lapply(1:length(df_names$name), function(i) 
  {
    print(i)
    plot_state_diagram(dataset_this=df_names$name[i], epsilon_this=e, capped_this = FALSE,
                       cost_add = 1,
                       cost_delete = 1,
                       cost_env = 1,
                       cost_wait = 0.1,
                       cex_vertex = 15.0 / df_names$`Number of species (n)`[i],
                       PLOT_SMALL = FALSE, PLOT_LABELS = FALSE)
    title(sprintf("(%s) %s - epsilon=%.1e", letters[i], df_names$Dataset[i], e),adj=0,cex.main=3,line=0)
  })
})
dev.off()






png(file='figures/g_state_diagram_invasion_all.png',width=10*3*100,height=10*2*100,res=100)
par(mfrow=c(2,3))
par(mar=c(0,0,4,0))
lapply(1:length(df_names$name), function(i) 
{
  print(i)
  plot_state_diagram(dataset_this=df_names$name[i], epsilon_this=0.1, capped_this = FALSE,
                     cost_add = 1,
                     cost_delete = 1,
                     cost_env = 1,
                     cost_wait = 0.1,
                     cex_vertex = 15.0 / df_names$`Number of species (n)`[i],
                     PLOT_SMALL = FALSE, PLOT_LABELS = FALSE,
                     DO_INVASION_GRAPH = TRUE)
  title(sprintf("(%s) %s", letters[i], df_names$Dataset[i]),adj=0,cex.main=3,line=0)
})
dev.off()




png(file='figures/g_state_diagram_uninvasion_all.png',width=10*3*100,height=10*2*100,res=100)
par(mfrow=c(2,3))
par(mar=c(0,0,4,0))
lapply(1:length(df_names$name), function(i) 
{
  print(i)
  plot_state_diagram(dataset_this=df_names$name[i], epsilon_this=0.1, capped_this = FALSE,
                     cost_add = 1,
                     cost_delete = 1,
                     cost_env = 1,
                     cost_wait = 0.1,
                     cex_vertex = 15.0 / df_names$`Number of species (n)`[i],
                     PLOT_SMALL = FALSE, PLOT_LABELS = FALSE,
                     DO_UNINVASION_GRAPH = TRUE)
  title(sprintf("(%s) %s", letters[i], df_names$Dataset[i]),adj=0,cex.main=3,line=0)
})
dev.off()





# make a big composite graph?
# TBD





### fig s1 distribution
plot_action_distribution <- function(dataset_this, 
                               epsilon_this, 
                               capped_this,
                               nice_name)
{
  df_aa_ss = df_aa %>% filter(epsilon==epsilon_this & dataset==dataset_this & capped==capped_this)
  
  fn_adj = df_aa_ss %>% filter(type=='adjacency') %>% pull(filename)
  
  transitions_this = read.csv(fn_adj)

  counts = transitions_this %>%
    mutate(temp = factor(temp, levels=c(0,1))) %>%
    group_by(add, del, temp) %>% 
    summarize(count=n())
  
  if (all(counts$temp==0))
  {
    counts = rbind(counts, data.frame(add=NA,del=NA,temp="1",count=NA))
  }
  
  ggplot(counts, aes(x=add,y=del,fill=log10(count))) + 
    geom_raster() + 
    coord_equal() + 
    scale_fill_viridis_c(name='Log10(Count)') +
    xlab('Number of additions') +
    ylab('Number of deletions') +
    theme_bw() +
    facet_wrap(~temp) +
    ggtitle(nice_name) +
    scale_y_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1)) +
    scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1))
}

plots_action_distribution = lapply(1:nrow(df_names), function(i) {
  plot_action_distribution(dataset_this=df_names$name[i], 
                           epsilon_this = 0.1, capped_this = FALSE,
                           nice_name=df_names$Dataset[i])
  })
g_plots_action_distribution = ggarrange(plotlist=plots_action_distribution, 
                                        common.legend = TRUE,legend = 'right',
                                        align='hv')
ggsave(g_plots_action_distribution, file='figures/g_plots_action_distribution.png',width=9,height=4)









# figure out most commonly visited states


df_all_by_group = df_all_joined %>% 
  filter(capped==FALSE) %>%
  rename(name=dataset) %>%
  left_join(df_names %>% select(name, Dataset)) %>%
  group_by(Dataset) %>% 
  group_split

counts_intermediate_visits = lapply(df_all_by_group, function(df) {
  print(df$Dataset[1])
  identities_top = rbindlist(lapply(1:nrow(df), function(i)
  {
    df_transitions_this = read.csv(df$filename_transitions[i])
    ids = df_transitions_this %>% 
      mutate(order=order(intermediate_visits, decreasing=TRUE)) %>%
      select(node,order) %>%
      pivot_wider(names_from=node, values_from=order)
    
    return(ids)
  }))
  
  return(identities_top)
})

g_counts_intermediate_visits_raw = lapply(1:length(counts_intermediate_visits), function(i)
{
  counts_long_this = counts_intermediate_visits[[i]] %>% 
    mutate(row=as.numeric(row.names(.))) %>% 
    cbind(df_all_by_group[[i]] %>% select(epsilon, cost_add, cost_delete, cost_environment, cost_wait)) %>%
    pivot_longer(!contains(c("row","epsilon","cost")),names_to='node',values_to='order') %>%
    mutate(log10_epsilon=paste("log10(epsilon)=",log10(epsilon)))
  
  g_counts_this = ggplot(counts_long_this,aes(x=row,y=node,fill=order<5)) +
    geom_raster() +
    scale_fill_manual(values=c('lightblue','red'),name='Most visited intermediate state (top 5)') +
    facet_wrap(~log10_epsilon,scales='free_x',ncol=5) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(df_all_by_group[[i]]$Dataset[1]) +
    xlab('A* experiment') +
    ylab('State')
})
g_counts_all_top_5_intermediate_visits_by_epsilon = ggarrange(plotlist=g_counts_intermediate_visits_raw, 
                                                              ncol=1,
                                                              common.legend = TRUE,
                                                              legend = 'bottom')
ggsave(g_counts_all_top_5_intermediate_visits_by_epsilon, file='figures/g_counts_all_top_5_intermediate_visits_by_epsilon.png',width=5,height=12)







# plot degree distribution
df_all_for_degree = df_all_joined %>% 
  rename(name=dataset) %>%
  left_join(df_names %>% select(name, Dataset)) %>%
  group_by(Dataset) %>% 
  filter(capped==FALSE, cost_add==1, cost_delete==1, cost_environment==1,cost_wait==1) %>%
  mutate(epsilon_pretty=paste("epsilon=",format(epsilon,scientific=FALSE)))
df_degree_all = rbindlist(lapply(1:nrow(df_all_for_degree), function(i) {
  z = read.csv(df_all_for_degree$filename_transitions[i]) %>%
    cbind(df_all_for_degree[i,,drop=FALSE])
  }))


g_degree = ggplot(df_degree_all, aes(x=in_degree,
                          y=out_degree, 
                          color=log10(intermediate_visits))) + 
  geom_point(alpha=0.75,size=0.5) + 
  theme_bw() + 
  scale_color_viridis_c(name="Log10 Number intermediate visits") +
  facet_grid(Dataset~epsilon_pretty) +
  scale_x_sqrt() +
  scale_y_sqrt() + 
  xlab("√In-degree") + 
  ylab("√Out-degree") +
  theme(legend.position='bottom')

ggsave(g_degree, file='figures/g_degree.png',width=8,height=9)














###### plot navigation probabilities
df_nav_probs = df_all_joined %>% 
  group_by(dataset,epsilon) %>% 
  filter(capped==FALSE) %>%
  summarize(nav_prob_mean = mean(Total.viable.path.proportions..viable.....candidate.pairs.), 
            nav_prob_sd=sd(Total.viable.path.proportions..viable.....candidate.pairs.)) %>%
  rename(name=dataset) %>%
  left_join(df_names %>% select(name, Dataset))

g_nav_prob = ggplot(df_nav_probs, aes(x=Dataset,
                          y=nav_prob_mean,
                          fill=factor(log10(epsilon)))) +
  geom_bar(stat='identity',position='dodge') +
  theme_bw() +
  ylim(0,1) +
  ylab('Navigation probability') +
  scale_fill_brewer(palette='Oranges',name=expression(paste("log"[10], "(", epsilon,")")))



# plot shortcut probability
# the index is Total.viable.paths..non.trivial.length. / Total.viable.paths
df_shortcut_probs = df_all_joined %>% 
  group_by(dataset,epsilon) %>% 
  filter(capped==FALSE) %>%
  summarize(shortcut_prob_mean = mean(Total.non.trivial.paths.proportions..non.trivial...viable.), 
            shortcut_prob_sd=sd(Total.non.trivial.paths.proportions..non.trivial...viable.)) %>%
  rename(name=dataset) %>%
  left_join(df_names %>% select(name, Dataset))

g_shortcut_prob = ggplot(df_shortcut_probs, aes(x=Dataset,
                                      y=shortcut_prob_mean,
                                      fill=factor(log10(epsilon)))) +
  geom_bar(stat='identity',position='dodge') +
  geom_errorbar(aes(ymin=shortcut_prob_mean-shortcut_prob_sd, 
                    ymax=shortcut_prob_mean+shortcut_prob_sd,
                    group=factor(log10(epsilon))),
                width=0.4,
                position=position_dodge(0.9)) +
  theme_bw() +
  ylim(0,1) +
  ylab('Shortcut probability') +
  scale_fill_brewer(palette='Purples',name=expression(paste("log"[10], "(", epsilon,")")))


# save combined figure
ggsave(ggarrange(g_nav_prob, g_shortcut_prob,
                 nrow=2,labels='auto',align='hv',
                 legend='right'), 
       file='figures/g_nav_shortcut_prob.png',width=9,height=7)





# assess path lengths
# make a histogram
df_length_counts = df_all_joined %>% 
  filter(capped==FALSE) %>%
  rename(name=dataset) %>%
  left_join(df_names %>% select(name, Dataset)) %>%
  select(contains("astar_length"),Dataset,epsilon) %>% 
  mutate(row=1:nrow(.)) %>%
  pivot_longer(!c("Dataset","row","epsilon"), names_to='length',values_to='fraction') %>% 
  mutate(length=as.numeric(gsub("astar_length_NOT_REACHABLE","",gsub("astar_length_X","",length)))) %>%
  filter(!is.na(length)) %>%
  filter(length > 1)

g_astar_length_counts = ggplot(df_length_counts, 
                               aes(x=as.numeric(as.character(factor(length))),
                                   y=fraction,group=row,color=factor(log10(epsilon)))) + 
  geom_line() + 
  geom_point(size=0.5) + 
  facet_wrap(~Dataset,scales='free_y') +
  theme_bw() +
  scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1),limits=c(2,8)) +
  scale_color_brewer(palette='Greens',name=expression(paste("log"[10], "(", epsilon,")"))) +
  theme(panel.grid.minor = element_blank()) +
  xlab("Shortcut path length") +
  ylab("Fraction") +
  scale_y_sqrt(limits=c(0,0.5))

ggsave(g_astar_length_counts, file='figures/g_astar_length_counts.png',width=7,height=4)









### PLOT PATHS

############################################################################################################
parse_names <- function(row_id)
{
  names_taxa_this = df_all_joined$Species.names[row_id]
  names_taxa_this = gsub("[","",names_taxa_this,fixed=TRUE)
  names_taxa_this = gsub("]","",names_taxa_this,fixed=TRUE)

  names_taxa_this = strsplit(names_taxa_this, split=", ")[[1]]
  names_taxa_this = sapply(names_taxa_this, function(x) { 
    x = sapply(strsplit(x, " ")[[1]], function(y) { substr(y, 1, 5) })
    paste(x, collapse=" ")
    })
  
  names_env_this = df_all_joined$Temperature.names[row_id]
  names_env_this = gsub("[","",names_env_this,fixed=TRUE)
  names_env_this = gsub("]","",names_env_this,fixed=TRUE)
  
  names_env_this = strsplit(names_env_this, split=", ")[[1]]
  if (length(names_env_this)==0)
  {
    names_env_this = ""
  }
  
  return(list(taxa=names_taxa_this,env=names_env_this))
}

replace_names <- function(path_string, parsed_names)
{
  if(substr(path_string,1,1) %in% c("+","-"))
  {
    action = substr(path_string,1,1)
    path_string = substr(path_string, 2, nchar(path_string)) # truncate
    path_string = strsplit(path_string, split=" +")[[1]]
    taxa = parsed_names$taxa[as.numeric(path_string)]
    final_string = paste(c(action, taxa),collapse="\n")
  }
  else if(substr(path_string,1,1) %in% c("*"))
  {
    action = substr(path_string,1,1)
    path_string = substr(path_string, 2, nchar(path_string)) # truncate

    env = parsed_names$env[as.numeric(path_string)]
    final_string = paste(c(action, env),collapse="\n")
  }
  else if(path_string==".")
  {
    final_string = path_string
  }
  else
  {
    path_string_split = strsplit(path_string, split="|", fixed=TRUE)[[1]]
    taxa = strsplit(path_string_split[1], split="*", fixed=TRUE)[[1]]
    env = strsplit(path_string_split[2], split="*", fixed=TRUE)[[1]]
    taxa_string = parsed_names$taxa[as.numeric(taxa)]
    env_string = parsed_names$env[as.numeric(env)]
    final_string = paste(c(taxa_string, env_string),collapse="\n")
    final_string = gsub("\n$","",final_string)
  }
  
  return(final_string)
}


plot_path <- function(path_this, parsed_names)
{
  if (is.na(path_this)) 
  {
    return(ggplot())
  }
  if (length(path_this)==0) 
  {
    return(ggplot())
  }
  # split out the actions and states
  path_parts_this = gsub("[\\{\\}]","",gsub(">","",strsplit(path_this, split='~~',fixed=TRUE)[[1]]))
  path_parts_this = unlist(strsplit(path_parts_this, ")(", fixed=TRUE))
  path_parts_this = gsub("(", "", path_parts_this, fixed=TRUE)
  # add in a wait action with a . symbol
  path_parts_this = gsub(")", ").", path_parts_this, fixed=TRUE)
  path_parts_this = unlist(strsplit(path_parts_this, ")", fixed=TRUE))
  
  # classify as states vs actions
  df_path_this = data.frame(path_part = path_parts_this)
  df_path_this$type = ifelse(substr(df_path_this$path_part,1,1)=="[","state","action")
  df_path_this$path_part = gsub("]", "", df_path_this$path_part,fixed=TRUE)
  df_path_this$path_part = gsub("[", "", df_path_this$path_part,fixed=TRUE)
  
  # specific classification
  df_path_this$action_type = ifelse(df_path_this$type=="state","state", 
                                    ifelse(substr(df_path_this$path_part,1,1)=='-',
                                           'deletion',
                                           ifelse(substr(df_path_this$path_part,1,1)=='+',
                                                  'addition',
                                                  ifelse(substr(df_path_this$path_part,1,1)=='.',
                                                         'wait',
                                                         ifelse(substr(df_path_this$path_part,1,1)=='*',
                                                                'environment',
                                                                'other')))))
  
  count_species_state <- function(path_part)
  {
    state = strsplit(path_part,split="|",fixed=TRUE)[[1]][1]
    state_parts = strsplit(state, "*", fixed=TRUE)[[1]]
    return(length(state_parts))
  }
  
  count_species_action <- function(path_part)
  {
    action_parts = strsplit(path_part, "[[:space:]]+")[[1]]
    return(length(action_parts))
  }
  
  
  df_path_this$richness = ifelse(df_path_this$type=="state", 
                                 sapply(df_path_this$path_part, count_species_state), NA)
  df_path_this$delta_richness = ifelse(df_path_this$type=="state", 
                                       0, ifelse(df_path_this$action_type=='addition', sapply(df_path_this$path_part,count_species_action),
                                                 ifelse(df_path_this$action_type=='deletion', -1*sapply(df_path_this$path_part,count_species_action),0)))
  
  df_path_this$sequence_count = NA
  df_path_this$delta_x = NA
  id_start = 1
  for (i in 1:nrow(df_path_this))
  {
    if (df_path_this$type[i]=='state')
    {
      df_path_this$sequence_count[i] = 0
      df_path_this$delta_x[id_start:i] = 1 / (1+max(df_path_this$sequence_count[id_start:i]))
      id_start = i+1
    }
    else
    {
      df_path_this$sequence_count[i] = df_path_this$sequence_count[i-1] + 1
    }
  }
  
  df_path_this$path_part_parsed = sapply(df_path_this$path_part, replace_names, parsed_names=parsed_names)
  
  #View(df_path_this)
  
  # add coordinates
  df_path_this = df_path_this %>% 
    fill(richness, .direction='down') %>%
    mutate(y = richness + delta_richness / max(delta_richness) / 4) %>%
    mutate(x = cumsum(delta_x) - 0.5) %>%
    mutate(action_type=factor(action_type, levels=c('state','addition','deletion','environment','wait'),ordered=TRUE))
  
  g = ggplot(df_path_this, aes(x=x,y=y,label=path_part_parsed,color=action_type)) +
    geom_line(data=df_path_this, aes(x=x,y=y), inherit.aes = FALSE, color='gray') +
    geom_label(aes(fill=action_type=='state'),size=2) + ### SIZE
    theme_bw() +
    scale_x_continuous(breaks=0:(max(df_path_this$x))) +
    scale_y_continuous(breaks=0:(max(df_path_this$richness)+1),limits=c(0,max(df_path_this$richness)+1)) +
    xlab("Action sequence") +
    ylab("Richness") +
    scale_color_manual(values=c("black","blue","red","orange","gray"),drop=FALSE) +
    scale_fill_manual(values=c('white','green')) +
    theme(legend.position='none') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(g)
}


# plot paths

plot_paths_example_indirect_removal = function(id_run, start_contains, goal_does_not_contain, which_id_ss, view=FALSE)
{
  paths_interesting = df_astar_paths_nontrivial[[id_run]] %>% 
    filter(a_star_ops > 1) %>% 
    # must lose a species but not directly
    filter(grepl(pattern=paste(start_contains), x=start, fixed=TRUE)==TRUE) %>% 
    filter(grepl(pattern=paste(goal_does_not_contain), x=goal, fixed=TRUE)==FALSE) %>% 
    filter(grepl(pattern=sprintf("(-%d)",goal_does_not_contain), x=full_path, fixed=TRUE)==FALSE) %>% 
    pull(full_path)
  
  print(length(paths_interesting))
  print(nrow(df_astar_paths_nontrivial[[id_run]]))
  
  if (view==TRUE)
  {
    View(data.frame(rowid=1:length(paths_interesting),paths_interesting))
  }
  
  metadata_interesting = df_all_joined %>% 
    slice(id_run) %>% 
    select(dataset:capped) # pick column header
  
  paths_interesting_ss = paths_interesting[which_id_ss]
  
  plots_paths_interesting_ss = lapply(1:length(paths_interesting_ss), function(i)
  {
    g = plot_path(paths_interesting_ss[i], parsed_names = parse_names(id_run))
  })
  
  return(list(g=plots_paths_interesting_ss,p=paths_interesting))
}


# find a clodif removal case
ids_bucci = df_all_joined %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(dataset=='Bucci' & 
           epsilon==0.1 & cost_add==1 & cost_delete==10 & cost_environment==1 & cost_wait==0.1 & capped==FALSE) %>%
  pull(rowid)

# 9 = clodif
g_paths_bucci_no_clodif = plot_paths_example_indirect_removal(ids_bucci, start_contains = 9, goal_does_not_contain=9, 
                                                              which_id_ss = c(92, 57))#c(92, 71, 57))



ids_maynard5 = df_all_joined %>% 
  mutate(rowid = 1:nrow(.)) %>%
  filter(dataset=='Maynard15-17-19-21-23' & 
           epsilon==0.1 & cost_add==1 & cost_delete==10 & cost_environment==1 & cost_wait==0.1 & capped==FALSE) %>%
  pull(rowid)




plot_paths_example_env_change = function(id_run, start_env, goal_env, which_id_ss=1, view=FALSE)
{
  paths_interesting = df_astar_paths_nontrivial[[id_run]] %>% 
    filter(a_star_ops > 1) %>% 
    filter(grepl(pattern="|1",x=start,fixed=TRUE)==TRUE) %>% # starting in environment one
    filter(grepl(pattern="|1",x=goal,fixed=TRUE)==TRUE) %>%
    filter(grepl(pattern="{(*",x=full_path,fixed=TRUE)==TRUE) %>%
    filter((net_add - net_del) < 0) %>%
    pull(full_path)
  
  print(length(paths_interesting))
  print(nrow(df_astar_paths_nontrivial[[id_run]]))
  
  if (view==TRUE)
  {
    View(data.frame(rowid=1:length(paths_interesting),paths_interesting))
  }
  
  metadata_interesting = df_all_joined %>% 
    slice(id_run) %>% 
    select(dataset:capped) # pick column header
  
  paths_interesting_ss = paths_interesting[which_id_ss]
  
  plots_paths_interesting_ss = lapply(1:length(paths_interesting_ss), function(i)
  {
    g = plot_path(paths_interesting_ss[i], parsed_names = parse_names(id_run))
  })

    return(list(g=plots_paths_interesting_ss,p=paths_interesting))
}



g_paths_maynard5 = plot_paths_example_env_change(ids_maynard5, which_id_ss = c(30,206))



g_paths = ggarrange(plotlist=c(g_paths_bucci_no_clodif$g, g_paths_maynard5$g), 
                    nrow=2,ncol=2,
                    labels='auto')
ggsave(g_paths, file='figures/g_paths.png', width=12,height=7)





















#### RANDOM FORESTS
calculate_jaccard <- function(from, to)
{
  # this works by ignoring the environmental variable (after the | in the string)
  substr_from = strsplit(from,"\\|")[[1]][1]
  substr_to = strsplit(to,"\\|")[[1]][1]
  
  vec_from = strsplit(substr_from, "\\*")[[1]]
  vec_to = strsplit(substr_to, "\\*")[[1]]
  
  index_jaccard = length(intersect(vec_from, vec_to)) / length(union(vec_from, vec_to))
  
  return(index_jaccard)
}

calculate_richness <- function(state)
{
  substr_state = strsplit(state,"\\|")[[1]][1]
  vec_state = strsplit(substr_state, "\\*")[[1]]
  
  return(length(vec_state))
}

# pick 100 paths from each parameter combination
df_all_transitions = rbindlist(lapply(1:nrow(df_all_joined), function(i) {
  cat(sprintf('%f\n',i/nrow(df_all_joined)))
  result = df_all_joined %>% 
    slice(i) %>% 
    select(dataset, capped, epsilon, starts_with("cost"), Proportion.of.candidate.states) %>% 
    cbind(df_astar_paths_nontrivial[[i]] %>% 
            select(-full_path) %>%
            sample_n(size=min(nrow(.), 100))) 
  
  return(result)
}))

### get a data subset ready for random forest modeling
# pick 20000 random slices per dataset (another downsample step)
df_all_transitions_ss = df_all_transitions %>% 
  group_by(dataset) %>%
  sample_n(min(nrow(df_all_transitions),20000)) %>%
  ungroup %>%
  mutate(path.exists = !is.na(a_star_ops)) %>%
  mutate(path.better.than.nominal = ifelse(path.exists, a_star_ops > 1, NA)) %>%
  mutate(path.type = ifelse(path.exists, ifelse(path.better.than.nominal, "shortcut", "direct"),"none"))

stats_transitions_ss = rbindlist(lapply(1:nrow(df_all_transitions_ss), function(i) {
  print(i/nrow(df_all_transitions_ss))
  
  jaccard = calculate_jaccard(from=df_all_transitions_ss$start[i], to=df_all_transitions_ss$goal[i])
  richness_from = calculate_richness(df_all_transitions_ss$start[i])
  richness_to = calculate_richness(df_all_transitions_ss$goal[i])
  richness_delta = richness_to - richness_from
  
  return(data.frame(jaccard=jaccard, 
                    richness_from = richness_from, 
                    richness_to=richness_to, 
                    richness_delta))
}))
df_all_transitions_ss = df_all_transitions_ss %>% 
  cbind(stats_transitions_ss)

# add dataset info for each dataset
df_one_each = df_all_joined %>% 
  group_by(dataset) %>% 
  slice(1) %>% 
  select(dataset, filename_summary)
df_one_each = df_one_each %>% cbind(rbindlist(lapply(df_one_each$filename_summary, function(fn) {
  df = read.csv(fn)[2:3,]
  df_split = strsplit(df, split=": ")
  return(data.frame(n=df_split[[1]][2], t=df_split[[2]][2]))
}))) %>%
  select(-filename_summary) %>%
  mutate(n=as.numeric(n), t=as.numeric(t))
df_all_transitions_ss = df_all_transitions_ss %>%
  left_join(df_one_each, by='dataset')

# add info on state
# assume the state info can be gotten from the first dataset/epsilon file (since this is just about the candidate states, not the dynamics
fns_assemblage = df_aa %>% filter(type=='assemblage') %>% 
  group_by(dataset, epsilon) %>% 
  slice(1) %>% 
  ungroup %>%
  select(filename, dataset, epsilon)
df_abundances_assemblage = rbindlist(lapply(1:nrow(fns_assemblage), function(i) {
  df = read.csv(fns_assemblage$filename[i]) %>%
    cbind(dataset=fns_assemblage$dataset[i], 
          epsilon=fns_assemblage$epsilon[i])
  return(df)
}))

df_all_transitions_ss = df_all_transitions_ss %>% mutate(primary_key_from = paste(dataset, epsilon, start), 
                                                         primary_key_to = paste(dataset, epsilon, goal)) %>%
  left_join(df_abundances_assemblage %>% 
              mutate(primary_key_from=paste(dataset, epsilon, str),
                     abundance_mean_from=abundance_mean) %>%
              select(primary_key_from, abundance_mean_from),
            by='primary_key_from') %>%
  left_join(df_abundances_assemblage %>% 
              mutate(primary_key_to=paste(dataset, epsilon, str),
                     abundance_mean_to=abundance_mean) %>%
              select(primary_key_to, abundance_mean_to),
            by='primary_key_to')

df_all_transitions_ss = df_all_transitions_ss %>%
  filter(capped==FALSE) %>%
  # NaNs are all where we have empty states
  mutate_at(vars(abundance_mean_from), ~replace(., is.nan(.), 0)) %>%
  mutate_at(vars(abundance_mean_to), ~replace(., is.nan(.), 0)) %>%
  # empty set jaccard similarities should be 1
  mutate_at(vars(jaccard), ~replace(., is.nan(.), 1)) %>%
  mutate(abundance_delta = abundance_mean_to - abundance_mean_from) %>%
  mutate(log10_epsilon = log10(epsilon)) %>%
  mutate(path.type = factor(path.type,levels=c('none','direct','shortcut')))

df_all_transitions_ss_balanced = df_all_transitions_ss %>%
  group_by(path.type) %>%
  sample_n(min(table(df_all_transitions_ss$path.type))) %>%
  ungroup


# make random forest models
m_rf_path_type = ranger(path.type ~ abundance_delta + log10_epsilon + richness_delta + jaccard + n + t + dataset + cost_add + cost_delete + cost_environment + cost_wait,
                        # balance the sample
                        data=df_all_transitions_ss_balanced,
                        importance='permutation',
                        probability=TRUE)

# check prediction error
m_rf_path_type$prediction.error

imp_raw = sort(importance(m_rf_path_type))
print(imp_raw)
df_imp = data.frame(var=names(imp_raw), value=imp_raw)
pred_names_nice = c(log10_epsilon="Log10(epsilon)",  # expression(paste("log"[10], "(", epsilon,")"))
                    n='Number of species (n)',
                    jaccard='Jaccard similarity',
                    dataset='Dataset',
                    t='Number of environments (m)',
                    abundance_delta='∆Abundance',
                    richness_delta='∆Richness',
                    cost_add='Cost (epsilon addition)',
                    cost_delete='Cost (epsilon deletion)',
                    cost_environment='Cost (environment)',
                    cost_wait='Cost (wait)')
df_imp$nice_name = pred_names_nice[df_imp$var]   

g_importance = ggplot(df_imp, aes(x=nice_name,y=value,fill=value)) + 
  geom_bar(stat='identity',fill=RColorBrewer::brewer.pal(3,'Blues')[3]) +
  theme_bw() +
  coord_flip() +
  ylab('Permutation importance') +
  xlab('Predictor') +
  theme(legend.position = 'none')
ggsave(g_importance, file='figures/g_importance.png',width=7,height=4)                    
                    
                    

# do pdps 
# this is slow - 5-10 minutes
pdps_path_type = rbindlist(lapply(1:length(levels(df_all_transitions_ss$path.type)), function(i) {
  print(i)
  df_pdp = pdp::partial(m_rf_path_type,
                        pred.var=c('richness_delta','abundance_delta','dataset'),
                        progress=TRUE,
                        type='classification',
                        which.class=i,
                        prob=TRUE,
                        grid.resolution=5)
  df_pdp = df_pdp %>% cbind(class=levels(df_all_transitions_ss$path.type)[i])
  
  return(df_pdp)
}))

g_rf_pdps = ggplot(pdps_path_type, aes(x=richness_delta,
                                          y=yhat,
                                          color=abundance_delta,
                                          group=abundance_delta)) + 
  geom_line(size=2,alpha=0.9) + 
  facet_grid(~class,scales='free_y') +
  scale_color_distiller(palette='Spectral',name='∆ Abundance',direction=1) +
  theme_bw() +
  xlab("∆ Richness") + 
  ylab("Probability of path type (predicted)") +
  ylim(0,1) +
  theme(legend.position='bottom')
ggsave(g_rf_pdps, file='figures/g_rf_pdps.png',width=7,height=5)



make_pdp_1d <- function(xvar)
{
  print(xvar)
  pdps_path_type = rbindlist(lapply(1:length(levels(df_all_transitions_ss$path.type)), function(i) {
    print(i)
    df_pdp = pdp::partial(m_rf_path_type,
                          pred.var=xvar,
                          progress=TRUE,
                          type='classification',
                          which.class=i,
                          prob=TRUE,
                          grid.resolution=5)
    df_pdp = df_pdp %>% cbind(class=levels(df_all_transitions_ss$path.type)[i])
    
    return(df_pdp)
  }))
  return(pdps_path_type)
}

pdps_1d_all = lapply(df_imp$var, make_pdp_1d)
# skip dataset as it is not continuous
pdps_1d_all_with_xvars = rbindlist(lapply(setdiff(1:length(pdps_1d_all),which(df_imp$var=="dataset")), function(i) {
  df_this = pdps_1d_all[[i]]
  df_this$xvar_name = pred_names_nice[ names(df_this)[1] ]
  names(df_this)[1] = "xvar"
  return(df_this)
  }))
  

g_pdps_1d = ggplot(pdps_1d_all_with_xvars, aes(x=xvar,y=yhat,color=class)) +
  facet_wrap(~xvar_name,scales='free',ncol=3) +
  geom_line(size=1, alpha=0.9) +
  theme_bw() +
  theme(legend.position='bottom') +
  ylim(0,1) +
  xlab("") + 
  ylab("Probability of path type (predicted)") +
  theme(axis.text.x=element_blank()) +
  scale_color_brewer(palette='Set1',name='Path type')
ggsave(g_pdps_1d, file='figures/g_pdps_1d.png', width=7,height=7)














## export table of species names
table_names = df_all_joined %>% 
  mutate(rowid=1:nrow(.)) %>%
  group_by(dataset, Species.names, Temperature.names) %>% 
  slice(1) %>%
  select(rowid, dataset, Species.names, Temperature.names)

table_names_full = rbindlist(lapply(1:length(table_names$rowid), function(i) { 
  taxa = parse_names(table_names$rowid[i])$taxa
  result = data.frame(name=table_names$dataset[i], number=1:length(taxa), abbreviation=taxa)
  return(result)
  })) %>%
  left_join(df_names %>% select(name, Dataset),by='name') %>%
  select(Dataset, everything()) %>%
  arrange(Dataset, number)

write.csv(table_names_full, file='figures/table_names_full.csv', row.names=FALSE)

table_env_full = rbindlist(lapply(1:length(table_names$rowid), function(i) { 
  env = parse_names(table_names$rowid[i])$env
  result = data.frame(name=table_names$dataset[i], number=1:length(env), abbreviation=env)
  return(result)
})) %>%
  left_join(df_names %>% select(name, Dataset),by='name') %>%
  select(Dataset, everything()) %>%
  arrange(Dataset, number)

write.csv(table_env_full, file='figures/table_env_full.csv', row.names=FALSE)






# check time constants
df_aa_assemblages = df_aa %>% filter(type=='assemblage')

taus = rbindlist(lapply(1:nrow(df_aa_assemblages), function(i) {
  df_assemblage_this = read.csv(df_aa_assemblages$filename[i])
  return(data.frame(tau=df_assemblage_this$tau,name=df_aa_assemblages$dataset[i]))
  })) %>%
  left_join(df_names %>% select(name, Dataset),by='name')

g_taus = ggplot(taus, aes(x=tau,fill=Dataset,color=Dataset)) + 
  geom_density(alpha=0.5) + 
  scale_x_log10() +
  theme_bw() +
  scale_color_brewer(palette='Set1') + 
  scale_fill_brewer(palette='Set1') +
  facet_wrap(~Dataset)

ggsave(g_taus,file='figures/g_taus.png',width=8,height=4)

