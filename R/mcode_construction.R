construction_mod <- function(input,local.net=FALSE,node.attribute=NULL,
                       db, species=c("human","ath"),
                       hierarchy=1)
{
  species <- match.arg(species)

  if(local.net){
    graph <- construct_local_mod(input=input,node.attribute=node.attribute)
  }else{
    species <- match.arg(species)
    graph <- construct_nlocal_mod(input=input,db=db,species=species,
                           hierarchy=hierarchy)
  }
  
  return(graph)
}


construct_local_mod <- function(input,node.attribute)
{
  if(missing(input)){
    stop("Not a data frame")
  }
  if(!is.data.frame(input)){
    stop("Not a data frame")
  }
  
  graph <- NULL
  
  if(nrow(input)>1){
    graph <- igraph::graph.data.frame(input,directed=FALSE)
    ##  add node attribute
    if(!is.null(node.attribute)){
      graph <- add.vertex.attri(graph=graph,input=node.attribute)
    }
  }else{
    stop("No edges are imported")
  }
  
  return(graph)
}


construct_nlocal_mod<-function(input,db, species=c("human","ath"),hierarchy=1)
{
  species <- match.arg(species)

  if(!missing(input)){
    if( !is.data.frame(input)){
      stop("Not a data frame")
    }
  }
  
     net <- igraph::graph.data.frame(db[,1:2],
                          directed=FALSE)
  
  
  if(missing(input)){
    return(net)
  }
  
  ##  match
  index <- match(input[,1],V(net)$name)
  index <- index[!is.na(index)]
  ##  create a  sub network
  graph <- igraph::induced.subgraph(graph=net, unlist(igraph::neighborhood(graph=net,order=hierarchy,nodes=index)))
  ##  vertex.hierarchy
  V(graph)$vertex.hierarchy<-rep(hierarchy,vcount(graph))
  if(hierarchy>=1){
    for(i in (hierarchy-1):0){
      index <- unlist(igraph::neighborhood(graph=graph,order=i,nodes=which(V(graph)$name %in% input[,1])))
      V(graph)$vertex.hierarchy[index]<-i
    }
  }
  ##  add vertex attributes
  graph <- add.vertex.attri(graph=graph,input=input)
  return(graph)
}

## Add attributes to the vertex. 
add.vertex.attri<-function(graph=NULL,input)
{
  if(!is.igraph(graph)){
    stop("Not an igraph object")
  }
  if(!is.data.frame(input)){
    input<-as.data.frame(input)
  }
  coln<-ncol(input)
  if(coln>1){
    colname<-colnames(input)
    index<-match(V(graph)$name,input[,1],nomatch=0)
    for(i in 2:coln){
      graph<-set.vertex.attribute(graph,colname[i],index=as.character(input[index,1]),
                                  value=input[index,i])
    }
  }
  return(graph)
}
