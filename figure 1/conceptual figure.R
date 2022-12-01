library(igraph)

set.seed(15)

el = as.matrix(data.frame(
  from=c("initial","initial", "t1","t2","t3", rep("initial",4)), 
  to=c("desired", "t1", "t2", "t3", "desired","o1","o2","o3","o4")))

g = graph_from_edgelist(el)
E(g)$type =c('direct',rep('shortcut',4),rep('brute-force',4))

l = layout_nicely(g)
l[2,1] = 0.5

pdf(file='g_illustration.pdf',width=5,height=5)
plot(g,
     layout=l,
     vertex.label=NA,
     edge.arrow.size=0.75,
     edge.color=ifelse(E(g)$type=='direct','black',ifelse(E(g)$type=='shortcut','gray','red')),
     vertex.color=ifelse(V(g)$name=='initial','purple',
                         ifelse(V(g)$name=='desired','green',
                                ifelse(grepl("^t",V(g)$name),'orange','lightgray'))),
     edge.color='black')
dev.off()