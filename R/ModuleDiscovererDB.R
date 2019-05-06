##< ##########################################################################################################################
##< # ModuleDiscoverer                                                                                                       #
##< # Copyright (c) 2015 Leibniz-Institut für Naturstoff-Forschung und Infektionsbiologie e.V. - Hans-Knoell-Institut (HKI). #
##< # Contributors: Sebastian Vlaic <Sebastian.Vlaic@hki-jena.de>                                                            #
##< # First version: November 2015                                                                                           #
##< # Filename: ModuleDiscovererDB.r                                                                                         #
##< # Language: R                                                                                                            #
##< #                                                                                                                        #
##< # This program is free software; you can redistribute it and/or modify it under the terms of the                         #
##< # GNU General Public License as published by the Free Software Foundation, Version 3.                                    #
##< #                                                                                                                        #
##< # This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied     #
##< # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  #
##< #                                                                                                                        #
##< #                                                                                                                        #
##< # ModuleDiscoverer is an algorithm for the identification of regulatory modules based on genome-wide, large-scale        #
##< # protein-protein interaction networks in conjunction with expression data from high-throughput experiments.             #
##< ##########################################################################################################################

moduleDiscoverer.fragmentGraph <- function(A=NULL, vlist=NULL, nbrOfSeeds=1, seed=NULL, verbose=FALSE){
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(A)){

    stop()
  }else{
    if(verbose){print(dim(A))}
  }
  
  if(is.null(vlist)){

    stop()
  }else{
    if(verbose){print(dim(vlist))}
  }
  
  if(!all(colnames(vlist) %in% c("weight","content","degree"))){
    stop()
  }
  
  if(length(unique(c(dim(A),dim(vlist)[1])))!=1){
    stop()
  }
  
  if(nbrOfSeeds>0){
  }else{
    stop()
  }
  
  # defines the vector of active nodes!
  currentNodes <- 1:dim(vlist)[1]
  
  # stores nodes that have to be excluded from A and vlist as final result
  exclude <- NULL
  
  # function to compute possible, valid neighbors
  # Only used for the initial search of minimal cliques of size 3.
  ## node: current seed node
  getNeighbors <- function(node){
    
    check <- intersect(which(A[node,]!=0), currentNodes)
    weights <- A[node,check]
    weight <- as.numeric(vlist[node,"weight"])
    
    check <- check[weights==weight]
    rm(weights)
    rm(weight)
    
    # re-order valid candidates!
    if(length(check)!=0){
      check <- check[sample(x=1:length(check), size=length(check))]
      # if the weight of the target node (the number of proteins it contains)
      # is unequal to the edge weight we remove the target node.
      # This happens if the target node is already a merged node!
      blacklist <- check[!(A[node,check]==as.numeric(vlist[check,"weight"]))]
      check <- setdiff(check, blacklist)
    }else{
      return(c())
    }
    
    return(check)
  }
  
  # identification of all three meres
  ## node: seed node
  find3mere <- function(node){
    # get all neighbors of the node
    check <- intersect(getNeighbors(node), currentNodes)
    
    if(length(check)>0){
      # for all neighbors check neighbors
      res = FALSE
      for(i in 1:length(check)){
        x = check[i]
        check2 <- setdiff(intersect(getNeighbors(x), currentNodes), node)
        if(length(check2)>0){
          for(j in 1:length(check2)){
            y = check2[j]
            tmp <- intersect(which(A[y,]!=0), currentNodes)
            if(any(tmp==node)){
              res <- TRUE
              names(res) <- paste(x,y,sep=".")
              break
            }
          }
        }
        if(res==TRUE){
          break;
        }
      }
      if(all(!res)){
        return(FALSE)
      }
    }else{
      return(FALSE)
    }
    return(as.numeric(strsplit(names(res)[1], split="\\.")[[1]]))
  }
  
  # extend the identified three meres
  ## node: contains one seed three meres
  extend3meres <- function(node){
    # get all neighbors of the clique
    check <- intersect(which(A[node,]!=0), currentNodes)
    
    # get the weights of the edges to all neighbors
    weights <- A[node,check]
    
    # get the weight of the node
    weight <- as.numeric(vlist[node,"weight"])
    
    # get the weights of all neighbors and multiply it by the weight of the node.
    # Why? This is necessary if two n-meres are merged. Consider two 3meres... if
    # two three meres are merged each node within one 3mere must have been initially
    # connected with all nodes of the other 3mere. Thus, the edge weight between the two
    # nodes must equal weight of the node times weight of the neighbor.
    nweight <- as.numeric(vlist[check,"weight"]) * weight
    
    # check only nodes for which their weight times the weight of the node equals the weight
    # of the edge that is connecting them!
    check <- check[weights==nweight]
    
    # If set by the user... randomize the order of the nodes to check!
    if(length(check)!=0){
      return(check[sample(x=1:length(check), size=length(check))][1])
    }else{
      return(FALSE)
    }
  }
  
  # merges at least two nodes into one! This will be the first node stored in nodes vector
  ## nodes: vector of nodes to be merged
  mergeNodes <- function(nodes){
    #for every node in nodes --> find neighbors and their weights!
    nbs <- lapply(nodes, function(node){
      nam <- intersect(which(A[node,]!=0), currentNodes)
      nbs <- A[node,nam]
      names(nbs) <- nam
      return(nbs)
    })
    if(verbose){print(paste("neighbors of",paste(nodes, collapse=", ")));print(nbs);print(paste("length of nodes vector:",length(nodes)))}
    
    
    # extract weights of neighbors
    weights <- rep(0, length(unique(names(unlist(nbs)))))
    names(weights) <- as.character(unique(names(unlist(nbs))))
    for(i in 1:length(nbs)){
      weights[names(nbs[[i]])] <- weights[names(nbs[[i]])] + as.numeric(nbs[[i]])
    }
    weights <- weights[as.character(setdiff(as.numeric(names(weights)),nodes))]
    names <- paste(vlist[nodes,"content"], collapse=" ")
    if(verbose){print(paste("names:",paste(names,collapse=", ")))}
    if(verbose){print(paste("weights:",paste(weights,collapse=", ")))}
    weight <- sum(as.numeric(vlist[nodes,"weight"]))
    
    # erase all but the first node from the nodes vector
    # this is actually a dangerous step since we alter the variables of the parent environment, but not of the current one!
    currentNodes <<- setdiff(currentNodes, nodes[-1]) # changes currentNodes in the parent environment!
    exclude <<- append(exclude, nodes[-1]) #this step is critical! don't do it anywhere else... # changes exclude in the parent environment!
    
    A[nodes[1],nodes[1]] <<- 0 # changes to A are made in the parent environment!
    A[nodes[1],as.numeric(names(weights))] <<- A[as.numeric(names(weights)),nodes[1]] <<- weights # changes to A are made to the parent environment!
    
    vlist[nodes[1],c("content","weight","degree")] <<- c(names, weight, sum(A[nodes[1],-append(exclude, nodes[-1])][A[nodes[1],-append(exclude, nodes[-1])]!=0])) # changes to A are made to the parent environment! Thats why we have to use the append(exclude, nodes[-1]) command again here. We have reset exclude earlier in the parent environment!
  }
  
  # this is were it all starts!
  
  if(verbose){print("computing all 3-meres...")}
  # step 1 merge all 3 meres... this takes O(n) worst case --> then the algorithm is done!
  tocheck <- sample(x=1:dim(vlist)[1], size=nbrOfSeeds)
  # we can actually exclude all nodes with out-degree 1 in the first search... (they will never form a 3-mere...)
  # also, this explains why the program always reports different numbers of nodes to check...
  tocheck <- tocheck[as.numeric(vlist[tocheck,"degree"])>1]
  
  if(verbose){print(paste("selected seed nodes:",paste(tocheck,collapse=", "),sep=" "))}
  
  if(verbose){print(paste("checking",length(tocheck),"nodes...",sep=" "))}
  if(length(tocheck)!=0){
    max <- length(tocheck)
    while(length(tocheck)!=0){
      if(verbose){print(paste("checking",tocheck[1]))}
      meres <- find3mere(tocheck[1])
      fuse <- NULL
      if(meres[1]){
        if(verbose){print(paste("fusing",tocheck[1],"with",paste(meres, collapse=", ")))}
        fuse <- c(tocheck[1],meres)
        mergeNodes(fuse)
      }else{
        if(verbose){print(paste("nothing found..."))}
        fuse <- tocheck[1]
      }
      tocheck <- setdiff(tocheck, fuse) # here we prevent that an identified 3-mere can be used as part of a new 3-mere ... This actually increases our chance to break out of very large cliques.
      if(verbose){print(paste("nodes to check:",paste(tocheck,collapse=", ")))}
      if(verbose){print("A matrix:");print(A[-exclude,-exclude]);print("vlist");print(vlist[-exclude,]);print("exclude");print(exclude);print("current nodes:");print(currentNodes)}
    }
  }
  if(verbose){print("done...")}
  
  
  
  if(verbose){print("extending 3 meres...")}
  if(length(exclude)==0){
    # now this can happen if no three-mere was formed in the first run. In this case we can stop here and return an empty list.
    return(vlist[-c(1:dim(vlist)[1]),,drop=FALSE])
  }
  # else we proceed...
  max <- dim(vlist)[1]-length(exclude)
  
  changed <- TRUE

  while(changed){
    changed <- FALSE
    # step 2 extend all 3 meres now worst case this can be O(n/3)
    tocheck <- sample(x=setdiff(1:dim(vlist)[1],exclude))
    
    # now we restrict the extension towards the identified k-meres. Why?
    # Well, because a single node can only be combined to another single node and build a 2-mere.
    # Otherwise it would have been identified in the 3-mere identification step.
    # However, starting from each 3-mere a single node can be added!
    # This logic has a crux if the vertices are initially partitioned, because this way 3-meres can
    # be missed in the first step. Given that the partition option should be used with random search
    # in a consensus approach there should be no problem... (Also we save CPU-runtime.... yea!)
    tocheck <- tocheck[vlist[tocheck,"weight"]>1]
    while(length(tocheck)>0){
      meres <- extend3meres(tocheck[1])
      if(meres[1]!=FALSE){
        fuse <- c(tocheck[1], meres)
        mergeNodes(fuse)
        changed <- TRUE
      }else{
        fuse <- tocheck[1]
      }
      tocheck <- setdiff(tocheck, fuse)
    }
   
  }
 
  
  vlist <- vlist[-exclude,,drop=FALSE]
  vlist <- vlist[as.numeric(vlist[,"weight"])>1,,drop=FALSE]
  return(vlist)
}

moduleDiscoverer.createDatabase <- function(results=NULL, proteins=NULL){
 
  if(is.null(results)){
    stop()
  }

  if(is.null(proteins)){
    stop()
  }

  proteinsInClique <<- c()
  results.indexed <- do.call(rbind, lapply(1:length(results), function(run){
    
    if(dim(results[[run]])[1]!=0){
      cbind("content"=sapply(results[[run]][,"content"], function(clique){
        x <- sort(as.numeric(unlist(strsplit(clique, split=" "))))
        proteinsInClique <<- unique(append(proteinsInClique, x))
        return(paste(x, collapse=" "))
      }), "run"=run)
    }else{
      return(NULL)
    }
  }))
  proteinsInClique <- unique(proteinsInClique)

  uniqueCliqueId <- unique(results.indexed[,1])
  names(uniqueCliqueId) <- uniqueCliqueId
  uniqueCliqueId[1:length(uniqueCliqueId)] = 1:length(uniqueCliqueId)
  results.indexed[,1] <- uniqueCliqueId[results.indexed[,1]]
  results.indexed <- apply(results.indexed,2,as.numeric)
  colnames(results.indexed)[1] <- "id"
  uniqueCliqueId <- names(uniqueCliqueId)

  database <- list("results"=results.indexed, "uniqueCliqueId"=uniqueCliqueId, "proteins"=proteins, "proteinsInClique"=proteinsInClique, "numberOfIterations"=length(results))

  return(database)
}

moduleDiscoverer.db.create_MD_object <- function(foregrounds=NULL, background=NULL, cores=1, chunks=5000, randomDataSets=NULL, minBackgroundRequired=3, minBackgroundToCliqueRatio=0.5, database=NULL, minForegroundRequired=1){
  if(is.null(database)){
     stop()
  }

  if(!is.null(foregrounds) & !is.list(foregrounds)){
     stop()
  }

  if(!is.numeric(cores) | cores < 1){
    stop()
  }

  if((!is.list(randomDataSets) & !is.null(randomDataSets)) | is.numeric(randomDataSets)){
    if(is.numeric(randomDataSets)){
      if(randomDataSets!=round(randomDataSets, digits=0)){
        stop()
      }else{
      }
    }else{

      stop()

      if(!all(unlist(lapply(randomDataSets, is.list)))){

        stop()
      }

      if(length(randomDataSets)!=length(foregrounds)){
       
        stop()
      }
    }
  }

  # components in the original interactome do not have to be necessarily
  # in a clique. Therefore, the background of the interactome is different
  # from the true background defined by the members of all cliques.
  components <- database$proteinsInClique
  names(components) <- database$proteins[components]

  if(is.null(background)){
    background <- components
  }else{
    tmp <- background
    background <- rep(-1, length(tmp))
    names(background) <- tmp
    rm(tmp)

    # and we filter the background according to the genes present in at least one clique...
    background[names(background) %in% names(components)] <- components[names(background)[names(background) %in% names(components)]]
  }

  background = background[background!=-1]

  fgs.names <- names(foregrounds)
  foregrounds = lapply(foregrounds, function(fgs){
    tmp = fgs
    fgs = rep(-1, length(tmp))
    names(fgs) = tmp
    rm(tmp)
    return(fgs)
  })

  numberOfForegrounds = length(foregrounds)
  if(!is.null(randomDataSets)){
    numberOfRandomDataSets <- unlist(lapply(randomDataSets, length))
  }else{
    numberOfRandomDataSets <- rep(0, numberOfForegrounds)
  }

  if(!is.null(randomDataSets) & !is.numeric(randomDataSets)){
    randomDataSets = unlist(randomDataSets, recursive=FALSE)
    randomDataSets = lapply(randomDataSets, function(rds){
      tmp = rds
      rds = rep(-1, length(tmp))
      names(rds) = tmp
      rm(tmp)
      return(rds)
    })
    foregrounds = append(foregrounds, randomDataSets)
  }

  foregrounds <- lapply(foregrounds, function(fgs){
    fgs[names(fgs) %in% names(components)] <- components[names(fgs)[names(fgs) %in% names(components)]]
    return(fgs)
  })

  foregrounds = lapply(foregrounds, function(fgs){
    return(fgs[fgs!=-1])
  })

  total <- length(database$uniqueCliqueId)

  if(is.numeric(randomDataSets)){
  
    randomDataSets <- lapply(1:numberOfForegrounds, function(i){
      return(lapply(1:randomDataSets, function(j){
        return(sample(x=background, size=length(foregrounds[[i]])))
      }))
    })
    numberOfRandomDataSets <- unlist(lapply(randomDataSets, length))
    randomDataSets <- unlist(randomDataSets, recursive=FALSE)
    foregrounds <- append(foregrounds, randomDataSets)

  }

 
  relevantCliques <- rep(-1, total)

  processCliques <- function(x){
    x <- as.numeric(unlist(strsplit(x, split=" ")))
    return(length(intersect(x,fgs))>=minForegroundRequired & length(intersect(x, background))>(minBackgroundRequired-1) & (length(intersect(x, background))/length(x))>minBackgroundToCliqueRatio)
  }

  fgs <- as.numeric(unique(unlist(sapply(1:numberOfForegrounds, function(i){return(foregrounds[[i]])}))))
  counter <- 1
  select <- rep(FALSE, total)
 
  while(counter<=total){
    end = counter+chunks-1
    if(end>total){
      end=total
    }
    
    runCores = cores
    if(parallel::detectCores()<=cores){
      runCores = parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(runCores, outfile = "")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    
    
    select[counter:end] <- foreach(clique = database$uniqueCliqueId[counter:end], .combine = 'c', .packages = "MODifieRDev") %dopar% {
      return(processCliques(clique))
    }
    counter <- counter+chunks

    if(counter==total){
      counter <- counter+1
    }
  }

  relevantCliques <- which(select)
  rm(fgs)
  relevantCliques <- setdiff(relevantCliques, -1)

  parallel::stopCluster(cl)

  if(!is.null(randomDataSets)){
    randomDataSets = lapply(setdiff(1:length(foregrounds),1:numberOfForegrounds), function(i){
      return(foregrounds[[i]])
    })

    randomDataSets = lapply(1:numberOfForegrounds, function(i){
      if(i==1){
        return(lapply(1:numberOfRandomDataSets[i], function(j){
          return(randomDataSets[[j]])
        }))
      }else{
        return(lapply((sum(numberOfRandomDataSets[1:(i-1)])+1):sum(numberOfRandomDataSets[1:i]), function(j){
          return(randomDataSets[[j]])
        }))
      }
    })
  }
  foregrounds = lapply(1:numberOfForegrounds, function(i){
    return(foregrounds[[i]])
  })
  names(foregrounds) <- fgs.names


  return(list("foregrounds"=foregrounds,"background"=background,"randomDataSets"=randomDataSets,"cores"=cores,"relevantCliques"=relevantCliques,"totalCliques"=total,"cores"=cores,"chunks"=chunks,"numberOfForegrounds"=numberOfForegrounds,"numberOfRandomDataSets"=numberOfRandomDataSets))
}

moduleDiscoverer.db.testForCliqueEnrichment <- function(database=NULL, input=NULL, cores=NULL, chunks=NULL){
  computeOnTheFly <- TRUE # Since there are usually a lot of cliques in combination with many random datasets setting this to FALSE is very likely of no use...

  if(is.null(database)){
    stop()
  }

  if(is.null(input)){
    stop()
  }

  if(is.null(cores)){
    cores <- input$cores
  }
  if(is.null(chunks)){
    chunks <- input$chunks
  }

  total <- input$totalCliques
  foregrounds <- input$foregrounds
  randomDataSets <- input$randomDataSets
  if(!is.null(randomDataSets)){
    foregrounds <- append(foregrounds, unlist(input$randomDataSets, recursive=FALSE))
  }
  background <- input$background
  relevantCliques <- input$relevantCliques
  numberOfForegrounds <- input$numberOfForegrounds
  numberOfRandomDataSets <- input$numberOfRandomDataSets

  runCores <- cores
  if(parallel::detectCores()<=cores){
      runCores <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(runCores)
  registerDoParallel(cl)
  
  parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  sets.n <- ceiling(length(relevantCliques)/chunks)
  sets <- list()
  for(i in 1:sets.n){
    tmp <- relevantCliques[((i-1)*chunks+1):(i*chunks)]
    tmp <- tmp[!is.na(tmp)]
    sets <- append(sets, list(tmp))
  }

  if(computeOnTheFly==TRUE){
    p.value <- matrix(NA, nrow=length(relevantCliques), ncol=numberOfForegrounds*2)
  }else{
    p.value <- matrix(NA, nrow=length(relevantCliques), ncol=length(foregrounds))
  }
  counter <- 0

  makeTest <- function(clique){
    members <- as.numeric(unlist(strsplit(clique, split=" ")))
    members <- intersect(members, background)
    res <- unlist(lapply(foregrounds, function(fgs){
      DC <- intersect(members, fgs)
      if(length(DC)!=0){
        DnC <- setdiff(fgs, members)
        nDC <- setdiff(members, DC)
        x <- matrix(c(length(DC),length(DnC),length(nDC),length(background)-length(DC)-length(DnC)-length(nDC)), byrow=T, ncol=2)
        m <- sum(x[, 1L])
        n <- sum(x[, 2L])
        k <- sum(x[1L, ])
        x <- x[1L, 1L]
        return(phyper(x - 1, m, n, k, lower.tail = FALSE))
      }else{
        return(1)
      }
    }))
    return(res)
  }

  makeTest.computeOnTheFly <- function(clique){
    members <- as.numeric(unlist(strsplit(clique, split=" ")))
    members <- intersect(members, background)
    res <- unlist(lapply(foregrounds, function(fgs){
      DC <- intersect(members, fgs)
      if(length(DC)!=0){
        DnC <- setdiff(fgs, members)
        nDC <- setdiff(members, DC)
        x <- matrix(c(length(DC),length(DnC),length(nDC),length(background)-length(DC)-length(DnC)-length(nDC)), byrow=T, ncol=2)
        m <- sum(x[, 1L])
        n <- sum(x[, 2L])
        k <- sum(x[1L, ])
        x <- x[1L, 1L]
        return(phyper(x - 1, m, n, k, lower.tail = FALSE))
      }else{
        return(1)
      }
    }))
    return(as.numeric(sapply(1:numberOfForegrounds, function(i){
      res.fgs <- res[i]
      res.rds <- res[-c(1:numberOfForegrounds)]
      if(numberOfRandomDataSets[i]!=0){
        if(i==1){
          res.rds <- res.rds[1:numberOfRandomDataSets[i]]
        }else{
          res.rds <- res.rds[(sum(numberOfRandomDataSets[1:(i-1)])+1):sum(numberOfRandomDataSets[1:i])]
        }
        return(c(sum(res.rds<=res.fgs)/length(res.rds),res.fgs))
      }else{
        return(c(NA, res.fgs))
      }
    })))
  }

  for(i in 1:length(sets)){
    ids <- sets[[i]]

    query.result <- cbind("id"=ids, "clique"=database$uniqueCliqueId[ids])

    if(computeOnTheFly==TRUE){
      p.value[(counter+1):(counter+dim(query.result)[1]),] <- foreach(i = 1:dim(query.result)[1], .combine = 'rbind') %dopar% {
        return(makeTest.computeOnTheFly(query.result[i,"clique"]))
      }
    }else{
      p.value[(counter+1):(counter+dim(query.result)[1]),] <- foreach(i = 1:dim(query.result)[1], .combine = 'rbind') %dopar% {
        return(makeTest(query.result[i,"clique"]))
      }
    }

    counter <- counter + dim(query.result)[1]
   
  }

  parallel::stopCluster(cl)
  if(!is.null(randomDataSets)){
    randomDataSets <- lapply(setdiff(1:length(foregrounds),1:numberOfForegrounds), function(i){
      return(foregrounds[[i]])
    })

    randomDataSets <- lapply(1:numberOfForegrounds, function(i){
      if(i==1){
        return(lapply(1:numberOfRandomDataSets[i], function(j){
          return(randomDataSets[[j]])
        }))
      }else{
        return(lapply((sum(numberOfRandomDataSets[1:(i-1)])+1):sum(numberOfRandomDataSets[1:i]), function(j){
          return(randomDataSets[[j]])
        }))
      }
    })
  }

  foregrounds <- lapply(1:numberOfForegrounds, function(i){
    return(foregrounds[[i]])
  })
  names(foregrounds) <- names(input$foregrounds)


  if(computeOnTheFly==TRUE){
    p.value <- lapply(seq(1,numberOfForegrounds*2,2), function(i){
      return(p.value[,c(i,i+1),drop=FALSE])
    })
  }else{
    pvals.fgs <- p.value[,1:numberOfForegrounds,drop=FALSE]
    pvals.rgs <- p.value[,setdiff(1:dim(p.value)[2],1:numberOfForegrounds),drop=FALSE]

    p.value <- lapply(1:numberOfForegrounds, function(i){
      if(i==1){
        return(cbind(pvals.fgs[,i], pvals.rgs[,1:numberOfRandomDataSets[i]]))
      }else{
        return(cbind(pvals.fgs[,i], pvals.rgs[,(sum(numberOfRandomDataSets[1:(i-1)])+1):sum(numberOfRandomDataSets[1:i])]))
      }
    })
  }

  return(list("foregrounds"=foregrounds, "background"=background, "randomDataSets"=randomDataSets, "p.value"=p.value, "cores"=cores, "numberOfForegrounds"=numberOfForegrounds, "numberOfRandomDataSets"=numberOfRandomDataSets, "relevantCliques"=relevantCliques, "computeOnTheFly"=computeOnTheFly, "chunks"=chunks))
}

moduleDiscoverer.db.extractEnrichedCliques <- function(database=NULL, result=NULL, p.value=0.05, coef=1, useFishersExactTestPvalue=FALSE){
  if(is.null(database)){
    stop()
  }
  if(is.null(result)){
  stop()
  }

  if(!(coef %in% 1:result$numberOfForegrounds)){
    stop()
  }

  if(all(is.na(result$p.value[[coef]][,1]))){
   
    ids <- as.numeric(result$relevantCliques[result$p.value[[coef]][,2]<p.value])
  }else{
    if(useFishersExactTestPvalue){
      ids <- as.numeric(result$relevantCliques[result$p.value[[coef]][,2]<p.value])
    }else{
      ids <- as.numeric(result$relevantCliques[result$p.value[[coef]][,1]<p.value])
    }
  }
  cliques <- list()
  if(length(ids)>0){
    query.result <- database$uniqueCliqueId[ids]

    cliques = append(cliques, lapply(1:length(query.result), function(i){
   
      members <- as.numeric(unlist(strsplit(query.result[i], split=" ")))
      members <- database$proteins[members]
      return(members)
    }))

   
    rm(query.result)
  }

  result$enrichedCliques <- cliques
  result$foregrounds <- result$foregrounds[[coef]]
  if(all(is.na(result$p.value[[coef]][,1]))){
    result$p.value <- result$p.value[[coef]][,2]
  }else{
    result$p.value <- result$p.value[[coef]][,1]
  }

  return(result)
}

moduleDiscoverer.module.createModule <- function(result=NULL, module.name="ModuleDiscoverer - regulatory module", cores=NULL){
  if(is.null(result)){
    stop()
  }

  if(is.null(cores)){
    cores <- result$cores
  }

  if(length(result$enrichedCliques)>0){

    runCores <- cores
    if(parallel::detectCores()<=cores){
     
      runCores <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(runCores)
    registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
    getEdges <- function(clique){
      return(t(combn(clique, m=2)))
    }
    ec <- result$enrichedCliques
    g <- igraph::simplify(graph.edgelist(el=foreach(i = 1:length(ec), .combine = 'rbind') %dopar% {
      return(getEdges(ec[[i]]))
    }, directed=FALSE))

    parallel::stopCluster(cl)

    g <- set.vertex.attribute(g, "inForeground", value=get.vertex.attribute(g)$name %in% names(result$foregrounds))
    g <- set.vertex.attribute(g, "inBackground", value=get.vertex.attribute(g)$name %in% names(result$background))
    tmp <- rep(0, length(get.vertex.attribute(g)$name))
    names(tmp) <- get.vertex.attribute(g)$name
    tmp[get.vertex.attribute(g)$name %in% names(result$background)] <- 1
    tmp[get.vertex.attribute(g)$name %in% names(result$foreground)] <- 2
    g <- set.vertex.attribute(g, "ForOrBacOrNon", value=tmp)
  }else{
    g <- NULL
  }

  return(g)
}
