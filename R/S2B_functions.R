##### All code in this file taken from S2B package and slightly modified
### https://github.com/frpinto/S2B

# This function allows you to get a vector of node ids matching your input 
# vector with strings of node identifiers. Identifiers not found in the 
# network are ignored.
seedrows=function(seed_graph,seedvec){
  vertexlist=unlist(igraph::vertex_attr(seed_graph))
  rowindex=which(is.element(vertexlist,seedvec))
}


simpmain=function(seed_graph){
  seed_graph_simp=igraph::simplify(seed_graph)
  seed_comp=igraph::components(seed_graph_simp)#divide el grafo en subgrafos
  main_id=which.max(seed_comp$csize)
  seed_graph_maincomp=igraph::induced_subgraph(seed_graph_simp,seed_comp[[1]]==main_id)
  seed_graph_maincomp
}

# S2B main function
# 
S2B_core = function(seed_graph,index1,index2,nrep,nrep2, n_cores = 4){
  meandist=igraph::mean_distance(seed_graph)
  bt=subS2B(seed_graph,index1,index2,meandist)
  pbt=rep(0,igraph::gorder(seed_graph))
  nscore=rep(0,igraph::gorder(seed_graph))
  deglist=igraph::degree(seed_graph)
  if (nrep2>0){
    rbt_matrix2=matrix(nrow=length(bt$allcount),ncol=nrep2)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    foreach(i = 1:nrep2,.packages = "MODifieRDev" ) %dopar% {
      
      rindex1=sample(igraph::gorder(seed_graph),length(index1),replace=FALSE)
      rindex2=sample(igraph::gorder(seed_graph),length(index2),replace=FALSE)
      rbt=subS2B(seed_graph,rindex1,rindex2,meandist)
      nscore[rbt$allcount<bt$allcount]=nscore[rbt$allcount<bt$allcount]+1
      rbt_matrix2[,i]=rbt$allcount
    }
    parallel::stopCluster(cl) # stop the cluster
    nscore=nscore/nrep2
  } else {
    rbt_matrix2=matrix()
  }
  if (nrep>0){
    rbt_matrix=matrix(nrow=length(bt$allcount),ncol=nrep)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    foreach(i = 1:nrep,.packages = "MODifieRDev" ) %dopar% {
      
      rg=igraph::sample_degseq(deglist,method="vl")
      rbt=subS2B(rg,index1,index2,meandist)
      pbt[rbt$allcount<bt$allcount]=pbt[rbt$allcount<bt$allcount]+1
      rbt_matrix[,i]=rbt$allcount
    }
    pbt=pbt/nrep
    parallel::stopCluster(cl) # stop the cluster
  } else {
    rbt_matrix=matrix()
  }
  bigvertexlist=igraph::vertex_attr(seed_graph)
  allstat=data.frame(protein=bigvertexlist[[1]],bcount=bt$allcount,score=pbt, nscore=nscore) ####
  s2btable=makes2btable(allstat,seed_graph,index1,index2)
  list(s2btable=s2btable,seedmat1=bt$smat1,seedmat2=bt$smat2,maxS2B=bt$maxS2B)
}

# S2B auxiliary function
#
# This function computes S2B scores, without specificity scores
subS2B=function(seed_graph,index1,index2,meandist){
  betweencount=rep(0,igraph::gorder(seed_graph)) #lista del tama??o de seed_graph
  seedmat1=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index1))
  seedmat2=matrix(data=0,nrow=igraph::gorder(seed_graph),ncol=length(index2))
  sp1=igraph::distances(seed_graph,v=index1,to=igraph::V(seed_graph))
  sp2=igraph::distances(seed_graph,v=index2,to=igraph::V(seed_graph))
  sp1[sp1==Inf]=igraph::vcount(seed_graph)
  sp2[sp2==Inf]=igraph::vcount(seed_graph)
  sp=sp1[,index2]
  maxbc=sum(sp>0 & sp<meandist)
  betweensub=union(index1,index2)
  for (i in 1:length(index1)){
    for (j in 1:length(index2)){
      m=sp1[i,index2[j]]
      if (m<meandist){
        sumsp=sp1[i,]+sp2[j,]
        nodelist=which(sumsp==m) # aquellos que cumplan la condicion, seran nodos presentes en un sh_path y seran a??adidos a nodelist
        betweencount[nodelist]=betweencount[nodelist]+1
        seedmat1[nodelist,i]=1
        seedmat2[nodelist,j]=1
      }else {
        nodelist=c(index1[i],index2[j])
        betweencount[nodelist]=betweencount[nodelist]+1
      }
    }
  }
  betweencount[index1]=betweencount[index1]-length(index2)
  betweencount[index2]=betweencount[index2]-length(index1)
  betweencount[intersect(index1,index2)]=betweencount[intersect(index1,index2)]+1
  betweencount=betweencount/maxbc
  list(allcount=betweencount,smat1=seedmat1,smat2=seedmat2,maxS2B=maxbc)
}

makes2btable=function(complete_results,apid_main,A_index,B_index){
  
  multi_neib2=multineibstat2(apid_main,A_index,B_index,1)
  #symblist=vertex_attr(Hs_uniprot2symb_graph(apid_main))
  s2b_apid_2=data.frame(id=complete_results[,1],S2B=complete_results[,2],S2Bspec=complete_results[,3],S2Bspec2=complete_results[,4],DA="candidate",D1neib=0,D2neib=0,bridges=0,bridgespec=0,stringsAsFactors=F)
  s2b_apid_2$DA[A_index]="Disease 1"
  s2b_apid_2$DA[B_index]="Disease 2"
  s2b_apid_2$DA[intersect(A_index,B_index)]="Disease 1/Disease 2"
  s2b_apid_2$D1neib[multi_neib2$neib]=multi_neib2$freq1
  s2b_apid_2$D2neib[multi_neib2$neib]=multi_neib2$freq2
  s2b_apid_2$bridges[multi_neib2$neib]=multi_neib2$prod
  s2b_apid_2$bridgespec[multi_neib2$neib]=multi_neib2$spec
  s2b_apid_2
  
}

multineibstat2=function(g,A_index,B_index,nrep){
  
  pn1=A_index
  for (i in 1:length(A_index)){
    pn1=c(pn1,igraph::neighbors(g,A_index[i]))
  }
  pn1_tab=tabulate(pn1,nbins=igraph::gorder(g))
  
  pn2=B_index
  for (i in 1:length(B_index)){
    pn2=c(pn2,igraph::neighbors(g,B_index[i]))
  }
  pn2_tab=tabulate(pn2,nbins=igraph::gorder(g))
  
  prod12=pn1_tab*pn2_tab
  
  pvec=rep(0,igraph::gorder(g))
  deglist=igraph::degree(g)
  for (k in 1:nrep){
    rg=igraph::sample_degseq(deglist,method="vl")
    rpn1=A_index
    for (i in 1:length(A_index)){
      rpn1=c(rpn1,igraph::neighbors(rg,A_index[i]))
    }
    rpn1_tab=tabulate(rpn1,nbins=igraph::gorder(g))
    
    rpn2=B_index
    for (i in 1:length(B_index)){
      rpn2=c(rpn2,igraph::neighbors(rg,B_index[i]))
    }
    rpn2_tab=tabulate(rpn2,nbins=igraph::gorder(g))
    
    rprod12=rpn1_tab*rpn2_tab
    pvec=pvec+(prod12>rprod12)+0.5*(prod12==rprod12)
  }
  pvec=pvec/nrep
  pn_multi=which(prod12>=1)
  data.frame(neib=pn_multi,freq1=pn1_tab[pn_multi],freq2=pn2_tab[pn_multi],spec=pvec[pn_multi],prod=prod12[pn_multi],stringsAsFactors=F)
}

s2bthreshold=function(s2bvec){
  n=length(s2bvec)
  x=max(s2bvec)*(1:n)/n
  y=sort(s2bvec)
  distvec=sqrt((x-max(s2bvec))^2+y^2)
  s2bt=y[which.min(distvec)]
  s2bt
}
