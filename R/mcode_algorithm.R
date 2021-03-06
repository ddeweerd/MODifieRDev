## vertex weighting
mcode.vertex.weighting<-function(graph,neighbors){	
  stopifnot(is.igraph(graph))  
  weight<-lapply(1:vcount(graph),
                 function(i){
                   subg<-induced.subgraph(graph,neighbors[[i]])
                   core<-graph.coreness(subg)
                   k<-max(core)
                   ### k-coreness
                   kcore<-induced.subgraph(subg,which(core==k))
                   if(vcount(kcore)>1){
                     if(any(is.loop(kcore))){
                       k*ecount(kcore)/choose(vcount(kcore)+1,2)						
                     }else{
                       k*ecount(kcore)/choose(vcount(kcore),2)
                     }
                   }else{
                     0
                   }
                 }
  )
  
  return(unlist(weight))
}

##	procedure MCODE-FIND-COMPLEX
mcode.find.complex<-function(neighbors,neighbors.indx,vertex.weight,
                             vwp,seed.vertex,seen)
{
  
  res <- .C("complex" ,as.integer(neighbors),as.integer(neighbors.indx),
          as.single(vertex.weight),as.single(vwp),as.integer(seed.vertex),
          seen=as.integer(seen),COMPLEX=as.integer(rep(0,length(seen)))
  , PACKAGE = "MODifieRDev")
  
  
  return(list(seen=res$seen,COMPLEX=which(res$COMPLEX!=0)))
}

## procedure MCODE-FIND-COMPLEXES
mcode.find.complexex<-function(graph,neighbors,vertex.weight,vwp)
{
  seen<-rep(0,vcount(graph))
  
  neighbors<-lapply(neighbors,function(item){item[-1]})
  neighbors.indx<-cumsum(unlist(lapply(neighbors,length)))
  
  neighbors.indx<-c(0,neighbors.indx)
  neighbors<-unlist(neighbors)-1
  
  COMPLEX<-list()
  n<-1
  w.order<-order(vertex.weight,decreasing=TRUE)
  for(i in w.order){
    if(!(seen[i])){
      res <- mcode.find.complex(neighbors,neighbors.indx,vertex.weight,vwp,i-1,seen)
      if(length(res$COMPLEX)>1){
        COMPLEX[[n]]<-res$COMPLEX
        seen<-res$seen
        n<-n+1
      }
    }
  }	
  rm(neighbors)	
  return(list(COMPLEX=COMPLEX,seen=seen))
}

## post-processing (optional)
mcode.fluff.complex<-function(graph,vertex.weight, fdt ,complex.g,seen)
{
  seq_complex.g<-seq_along(complex.g)
  for(i in seq_complex.g){
    node.neighbor<-unlist(neighborhood(graph,1,complex.g[i]))
    if(length(node.neighbor)>1){
      subg<-induced.subgraph(graph,node.neighbor)
      if(graph.density(subg, loops=FALSE)>fdt){
        complex.g<-c(complex.g,node.neighbor)
      }
    }
  }
  
  return(unique(complex.g))
}

#	procedure MCODE-POST-PROCESS
mcode.post.process<-function(graph,vertex.weight,haircut,fluff, fdt,
                             set.complex.g,seen)
{
  indx<-unlist(lapply(set.complex.g,
                      function(complex.g){
                        if(length(complex.g)<=2)
                          0
                        else
                          1
                      }
  ))
  set.complex.g<-set.complex.g[indx!=0]
  set.complex.g<-lapply(set.complex.g,
                        function(complex.g){
                          coreness<-graph.coreness(induced.subgraph(graph,complex.g))						
                          if(fluff){
                            complex.g<-mcode.fluff.complex(graph,vertex.weight,fdt,complex.g,seen)
                            if(haircut){
                              ## coreness needs to be recalculated
                              coreness<-graph.coreness(induced.subgraph(graph,complex.g))
                              complex.g<-complex.g[coreness>1]
                            }
                          }else if(haircut){
                            complex.g<-complex.g[coreness>1]
                          }
                          return(complex.g)
                        })
  set.complex.g<-set.complex.g[lapply(set.complex.g,length)>2]
  return(set.complex.g)
}

mcode <- function(graph, vwp=0.5, haircut=FALSE, fluff=FALSE, fdt, loops=TRUE)
{
  stopifnot(is.igraph(graph))

  if(!loops){
    graph<-simplify(graph,remove.multiple=FALSE,remove.loops=TRUE)
  }
  neighbors<-neighborhood(graph,1)
  W<-mcode.vertex.weighting(graph,neighbors)
  res<-mcode.find.complexex(graph,neighbors=neighbors,vertex.weight=W,vwp=vwp)
  COMPLEX<-mcode.post.process(graph,vertex.weight=W,haircut=haircut,fluff=fluff,
                              fdt=fdt,res$COMPLEX,res$seen)		
  score<-unlist(lapply(COMPLEX,
                       function(complex.g){
                         complex.g<-induced.subgraph(graph,complex.g)
                         if(any(is.loop(complex.g)))
                           score<-ecount(complex.g)/choose(vcount(complex.g)+1,2)*vcount(complex.g)
                         else
                           score<-ecount(complex.g)/choose(vcount(complex.g),2)*vcount(complex.g)
                         return(score)
                       }
  ))
  order_score<-order(score,decreasing=TRUE)
  return(list(COMPLEX=COMPLEX[order_score],score=score[order_score]))
}