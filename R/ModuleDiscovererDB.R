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

#' @export
moduleDiscoverer.fragmentGraph <- function(A=NULL, vlist=NULL, nbrOfSeeds=1, seed=NULL, verbose=FALSE){
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(A)){
    cat(paste('Adjacency matrix has to be given!','\n'))
    stop()
  }else{
    if(verbose){print(dim(A))}
  }
  
  if(is.null(vlist)){
    cat(paste('Vertex list has to be given!','\n'))
    stop()
  }else{
    if(verbose){print(dim(vlist))}
  }
  
  if(!all(colnames(vlist) %in% c("weight","content","degree"))){
    cat(paste('Cols of vertex list have to include "weight", "content", "degree"!','\n'))
    stop()
  }
  
  if(length(unique(c(dim(A),dim(vlist)[1])))!=1){
    cat(paste("Vertex matrix's nrow must equal dim of adjacency matrix!",'\n'))
    stop()
  }
  
  if(nbrOfSeeds>0){
  }else{
    cat(paste('nbrOfSeeds has to be > 0!','\n'))
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
  if(verbose){pb <- txtProgressBar(min=0, max=max, initial=0, style=3)}
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
    if(verbose){setTxtProgressBar(pb, value=max-length(tocheck))}
  }
  if(verbose){close(pb)}
  
  vlist <- vlist[-exclude,,drop=FALSE]
  vlist <- vlist[as.numeric(vlist[,"weight"])>1,,drop=FALSE]
  return(vlist)
}





moduleDiscoverer.createDatabase <- function(results=NULL, proteins=NULL){
 
  if(is.null(results)){
    cat(paste('results is null','\n'))
    stop()
  }

  if(is.null(proteins)){
    cat(paste('proteins is null','\n'))
    stop()
  }

  cat(paste('INFO: indexing results...','\n'))
  #pb <- txtProgressBar(min = 0, max = length(results), style = 3)
  proteinsInClique <<- c()
  results.indexed <- do.call(rbind, lapply(1:length(results), function(run){
    #setTxtProgressBar(pb, run)
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
  #close(pb)

  cat(paste('INFO: indexing cliques...','\n'))
  uniqueCliqueId <- unique(results.indexed[,1])
  names(uniqueCliqueId) <- uniqueCliqueId
  uniqueCliqueId[1:length(uniqueCliqueId)] = 1:length(uniqueCliqueId)
  results.indexed[,1] <- uniqueCliqueId[results.indexed[,1]]
  results.indexed <- apply(results.indexed,2,as.numeric)
  colnames(results.indexed)[1] <- "id"
  uniqueCliqueId <- names(uniqueCliqueId)

  cat(paste('INFO: preparing results...','\n'))
  database <- list("results"=results.indexed, "uniqueCliqueId"=uniqueCliqueId, "proteins"=proteins, "proteinsInClique"=proteinsInClique, "numberOfIterations"=length(results))

  cat(paste('INFO: all done!','\n'))
  return(database)
}

moduleDiscoverer.db.create_MD_object <- function(foregrounds=NULL, background=NULL, cores=1, chunks=5000, randomDataSets=NULL, minBackgroundRequired=3, minBackgroundToCliqueRatio=0.5, database=NULL, minForegroundRequired=1){
  if(is.null(database)){
    cat(paste('ERROR: no database supplied!','\n'))
    stop()
  }

  if(!is.null(foregrounds) & !is.list(foregrounds)){
    cat(paste('ERROR: no list of foregrounds provided','\n'))
    stop()
  }

  if(!is.numeric(cores) | cores < 1){
    cat(paste("ERROR: cores has to be a positive integer value!",'\n'))
    stop()
  }

  if((!is.list(randomDataSets) & !is.null(randomDataSets)) | is.numeric(randomDataSets)){
    if(is.numeric(randomDataSets)){
      if(randomDataSets!=round(randomDataSets, digits=0)){
        cat(paste("ERROR: randomDataSets has to be an integer > 0",'\n'))
        stop()
      }else{
        cat(paste("INFO:",randomDataSets,"random data sets will be created",'\n'))
      }
    }else{
      cat(paste("ERROR: randomDataSets has to be a list of lists of random gene sets",'\n'))
      stop()

      if(!all(unlist(lapply(randomDataSets, is.list)))){
        cat(paste("ERROR: all elements of randomDataSets have to be a list of random data sets",'\n'))
        stop()
      }

      if(length(randomDataSets)!=length(foregrounds)){
        cat(paste("ERROR: the length of randomDataSets has to equal the number of passed foregrounds",'\n'))
        stop()
      }
    }
  }

  # components in the original interactome do not have to be necessarily
  # in a clique. Therefore, the background of the interactome is different
  # from the true background defined by the members of all cliques.
  cat(paste("INFO: reducing the interactome background...",'\n'))
  components <- database$proteinsInClique
  names(components) <- database$proteins[components]
  cat(paste('\t',"INFO: there are",length(database$proteins),"ids in the interactome.",'\n'))
  cat(paste('\t',"INFO:",length(components),"ids are part of a clique.",'\n'))

  cat(paste("INFO: mapping background to internal IDs...",'\n'))
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

  cat(paste('\t',"INFO: done, ",length(background)-sum(background==-1)," of ",length(background)," components were successfully mapped",'\n',sep=""))
  background = background[background!=-1]

  cat(paste("INFO: mapping foreground to internal IDs... "))
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

  cat(paste("done!",'\n',sep=""))
  for(i in 1:numberOfForegrounds){
    if(!is.null(names(foregrounds))){
      cat(paste('\t',"INFO: Foreground ",i," (",names(foregrounds)[[i]],"): ",length(foregrounds[[i]])-sum(foregrounds[[i]]==-1)," of ",length(foregrounds[[i]])," components were successfully mapped","\n",sep=""))
    }else{
      cat(paste('\t',"INFO: Foreground ",i,": ",length(foregrounds[[i]])-sum(foregrounds[[i]]==-1)," of ",length(foregrounds[[i]])," components were successfully mapped","\n",sep=""))
    }
  }

  foregrounds = lapply(foregrounds, function(fgs){
    return(fgs[fgs!=-1])
  })

  total <- length(database$uniqueCliqueId)

  if(is.numeric(randomDataSets)){
    cat(paste("INFO: creating random data sets..."))
    randomDataSets <- lapply(1:numberOfForegrounds, function(i){
      return(lapply(1:randomDataSets, function(j){
        return(sample(x=background, size=length(foregrounds[[i]])))
      }))
    })
    numberOfRandomDataSets <- unlist(lapply(randomDataSets, length))
    randomDataSets <- unlist(randomDataSets, recursive=FALSE)
    foregrounds <- append(foregrounds, randomDataSets)
    cat(paste('done','\n'))
  }

  cat(paste("INFO: initializing parallel environment...",'\n'))
  runCores = cores
  if(detectCores()<=cores){
    cat(paste('\n','\t',"WARNING: Maximal number of cores detected on the system is ",detectCores(),". I'll use ",detectCores()-1," cores instead!",'\n',sep=""))
    runCores = detectCores() - 1
  }
  cl <- makeCluster(runCores)
  registerDoParallel(cl)


  cat(paste("INFO: identifying relevant cliques... ",'\n'))
  relevantCliques <- rep(-1, total)

  processCliques <- function(x){
    x <- as.numeric(unlist(strsplit(x, split=" ")))
    return(length(intersect(x,fgs))>=minForegroundRequired & length(intersect(x, background))>(minBackgroundRequired-1) & (length(intersect(x, background))/length(x))>minBackgroundToCliqueRatio)
  }

  fgs <- as.numeric(unique(unlist(sapply(1:numberOfForegrounds, function(i){return(foregrounds[[i]])}))))
  counter <- 1
  select <- rep(FALSE, total)
  #pb <- txtProgressBar(min = counter, max = total, style = 3)
  while(counter<=total){
    end = counter+chunks-1
    if(end>total){
      end=total
    }
    select[counter:end] <- foreach(clique = database$uniqueCliqueId[counter:end], .combine = 'c') %dopar% {
      return(processCliques(clique))
    }
    counter <- counter+chunks
    #setTxtProgressBar(pb, counter)
    if(counter==total){
      counter <- counter+1
    }
  }
  #close(pb)
  relevantCliques <- which(select)
  rm(fgs)
  relevantCliques <- setdiff(relevantCliques, -1)

  cat(paste("INFO: stopping parallel environment...",'\n'))
  stopCluster(cl)

  cat(paste("INFO: distangling data...",'\n'))
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

  cat(paste("INFO: all done!",'\n'))

  return(list("foregrounds"=foregrounds,"background"=background,"randomDataSets"=randomDataSets,"cores"=cores,"relevantCliques"=relevantCliques,"totalCliques"=total,"cores"=cores,"chunks"=chunks,"numberOfForegrounds"=numberOfForegrounds,"numberOfRandomDataSets"=numberOfRandomDataSets))
}

moduleDiscoverer.db.testForCliqueEnrichment <- function(database=NULL, input=NULL, cores=NULL, chunks=NULL){
  computeOnTheFly <- TRUE # Since there are usually a lot of cliques in combination with many random datasets setting this to FALSE is very likely of no use...

  if(is.null(database)){
    cat(paste('ERROR: no database supplied!','\n'))
    stop()
  }

  if(is.null(input)){
    cat(paste('ERROR: input is null','\n'))
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

  cat(paste("INFO: initializing parallel environment...",'\n'))
  runCores <- cores
  if(detectCores()<=cores){
    cat(paste('\t',"WARNING: Maximal number of cores detected on the system is ",detectCores(),". Using ",detectCores()-1," cores instead!",'\n',sep=""))
    runCores <- detectCores() - 1
  }
  cl <- makeCluster(runCores)
  registerDoParallel(cl)
  cat(paste("INFO: done!",'\n'))

  sets.n <- ceiling(length(relevantCliques)/chunks)
  sets <- list()
  for(i in 1:sets.n){
    tmp <- relevantCliques[((i-1)*chunks+1):(i*chunks)]
    tmp <- tmp[!is.na(tmp)]
    sets <- append(sets, list(tmp))
  }

  cat(paste("INFO: testing ",length(relevantCliques)," (of ",total," in total) relevant cliques for enrichment...",'\n',sep=""))
  cat(paste("INFO: initializing variables...",'\n'))
  if(computeOnTheFly==TRUE){
    cat(paste('\t',"INFO: using permutation based p-values calculation on the fly...",'\n'))
    p.value <- matrix(NA, nrow=length(relevantCliques), ncol=numberOfForegrounds*2)
  }else{
    p.value <- matrix(NA, nrow=length(relevantCliques), ncol=length(foregrounds))
  }
  cat(paste("INFO: done!",'\n'))
  #pb <- txtProgressBar(min = 0, max = length(relevantCliques), style = 3)
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
    #setTxtProgressBar(pb, counter)
  }
  #close(pb)
  #rm(pb)
  cat(paste("INFO: done!",'\n'))

  cat("INFO: stopping parallel environment...")
  stopCluster(cl)
  cat(paste("done!",'\n'))

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
    cat(paste('ERROR: no database supplied!','\n'))
    stop()
  }
  if(is.null(result)){
    cat(paste('ERROR: result is null','\n'))
    stop()
  }

  if(!(coef %in% 1:result$numberOfForegrounds)){
    cat(paste('ERROR: coef has to be between 1 and ',result$numberOfForegrounds,'\n'))
    stop()
  }

  cat(paste('INFO: extracting significantly enriched cliques...','\n'))
  if(all(is.na(result$p.value[[coef]][,1]))){
    cat(paste('\t','WARNING: no random datasets were supplied. Using Fisher\'s exact test p-values instead...','\n'))
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

    #pb <- txtProgressBar(min = 1, max = length(query.result), style = 3)
    cliques = append(cliques, lapply(1:length(query.result), function(i){
      #setTxtProgressBar(pb, i)
      members <- as.numeric(unlist(strsplit(query.result[i], split=" ")))
      members <- database$proteins[members]
      return(members)
    }))

    #close(pb)
    rm(query.result)
  }

  result$enrichedCliques <- cliques
  result$foregrounds <- result$foregrounds[[coef]]
  if(all(is.na(result$p.value[[coef]][,1]))){
    result$p.value <- result$p.value[[coef]][,2]
  }else{
    result$p.value <- result$p.value[[coef]][,1]
  }

  cat(paste('INFO: all done!','\n'))
  return(result)
}

moduleDiscoverer.db.testNetworkEnrichment <- function(database=NULL, result=NULL, p.value=0.05, chunks=NULL, cores=NULL){

  if(is.null(database)){
    cat(paste('ERROR: no database supplied!','\n'))
    stop()
  }

  if(is.null(result)){
    cat(paste('ERROR: result is null','\n'))
    stop()
  }

  if(is.null(result$randomDataSets)){
    cat(paste('WARNING: no random dataset supplied','\n'))
    pvalColumn <- 2
  }else{
    pvalColumn <- 1
  }

  if(is.null(cores)){
    cores <- result$cores
  }
  if(is.null(chunks)){
    chunks <- result$chunks
  }

  cat(paste("INFO: initializing parallel environment...",'\n'))
  runCores <- cores
  if(detectCores()<=cores){
    cat(paste('\t',"WARNING: Maximal number of cores detected on the system is ",detectCores(),". Using ",detectCores()-1," cores instead!",'\n',sep=""))
    runCores <- detectCores() - 1
  }
  cl <- makeCluster(runCores)
  registerDoParallel(cl)
  cat(paste("INFO: done!",'\n'))

  networkEnrichment = matrix(NA, ncol=result$numberOfForegrounds, nrow=2)
  rownames(networkEnrichment) = c("random datasets (statistical background)","Fisher's exact test (whole interactome)")
  if(all(!is.na(names(result$foregrounds)))){
    colnames(networkEnrichment) = names(result$foregrounds)
  }else{
    colnames(networkEnrichment) = paste("foreground-",1:result$numberOfForegrounds,sep="")
  }
  for(j in 1:result$numberOfForegrounds){
    if(is.null(names(result$foregrounds))){
      cat(paste('INFO: processing foreground',j,'\n'))
    }else{
      cat(paste('INFO: processing foreground (',names(result$foregrounds)[j],')',j,'\n'))
    }

    cat(paste('INFO: creating full-set result...'))
    ids <- result$relevantCliques[result$p.value[[j]][,pvalColumn]<p.value]
    query.result <- database$uniqueCliqueId[ids]
    nodes <- unique(unlist(lapply(query.result, function(x){unlist(strsplit(x, split=" "))})))
    cat(paste('done!','\n'))

    DM <- sum(nodes %in% result$foregrounds[[j]]) # number of DEGs in the the regulatory module
    nDM <- sum(!(nodes %in% result$foregrounds[[j]])) # number of non DEGs in the regulatory module. This includes also proteins that are not in the statistical background of the high-throughput platform!!!
    DnM <- sum(!(result$foregrounds[[j]] %in% nodes)) # number of DEGs that are not in the regulatory module
    nDnM <- length(database$proteins[database$proteinsInClique]) - DM - nDM - DnM
    ftpval <- fisher.test(matrix(c(DM,nDM,DnM,nDnM), ncol=2, byrow=T), alternative="greater")$p.value

    if(result$numberOfRandomDataSets[j]!=0){
      cat(paste('INFO: processing randomDataSets...','\n'))
      #pb <- txtProgressBar(min = 1, max = result$numberOfRandomDataSets[j], style = 3)
      tmp <- rep(-1, result$numberOfRandomDataSets[j])
      counter <- 1
      while(counter<=result$numberOfRandomDataSets[j]){
        end <- counter+chunks-1
        if(end>result$numberOfRandomDataSets[j]){
          end <- result$numberOfRandomDataSets[j]
        }
        tmp[counter:end] <- foreach(i = counter:end, .combine = 'c') %dopar% {
          return(sum(result$randomDataSets[[j]][[i]] %in% nodes))
        }
        counter <- counter+chunks
        #setTxtProgressBar(pb, counter)
        if(counter==result$numberOfRandomDataSets[j]){
          counter <- counter + 1
        }
      }
      #close(pb)
      cat(paste('done!','\n'))
      tmp <- 1-(sum(tmp<=(sum(result$foregrounds[[j]] %in% nodes)))/result$numberOfRandomDataSets[j])
      if(tmp==0){
        tmp = 1/result$numberOfRandomDataSets[j]
      }

      networkEnrichment[,j] <- c(tmp,ftpval)
    }else{
      cat(paste("WARNING: no random datasets are supplied for this foreground",'\n'))
      networkEnrichment[2,j] <- ftpval
    }
  }

  cat("INFO: stopping parallel environment...")
  stopCluster(cl)
  cat(paste("done!",'\n'))

  return(networkEnrichment)
}

moduleDiscoverer.db.computeNetworkStability <- function(database=NULL, result=NULL, repeats=100, p.value=0.05, chunks=NULL, cores=NULL, size=NULL, module.reference=NULL){

  if(is.null(database)){
    cat(paste('ERROR: no database supplied!','\n'))
    stop()
  }

  total <- database$numberOfIterations

  if(is.null(result)){
    cat(paste('ERROR: result is null','\n'))
    stop()
  }

  if(is.null(size)){
    size <- total
  }else if(size>total){
    cat(paste("INFO: Size is larger than total number iterations. Using total number of iterations instead!",'\n'))
    size <- total
  }

  if(is.null(cores)){
    cores <- result$cores
  }
  if(is.null(chunks)){
    chunks <- result$chunks
  }

  if(!is.null(module.reference)){
    cat(paste("INFO: a reference graph was provided, calculating difference to reference solution...",'\n'))
    if(class(module.reference)!="igraph"){
      cat(paste("ERROR: reference module has to be an igraph object.",'\n'))
      stop()
    }
  }

  cat(paste("INFO: initializing parallel environment...",'\n'))
  runCores <- cores
  if(detectCores()<=cores){
    cat(paste('\t',"WARNING: Maximal number of cores detected on the system is ",detectCores(),". I'am initializing accordingly (max cores - 1)!",'\n',sep=""))
    runCores <- detectCores() - 1
  }
  cl <- makeCluster(runCores)
  registerDoParallel(cl)
  cat(paste("INFO: done!",'\n'))

  cat(paste('INFO: creating',repeats,'sub-sets for',size,'iterations...','\n'))
  bootstrap.samples <- lapply(1:repeats, function(k){
    return(sample(unique(database$results[,"run"]), replace=TRUE, size=size))
  })

  subsampling.results  <- matrix(NA, ncol=8, nrow=0)

  rep.nodes <- list()
  rep.edges <- list()
  rep.fuSet.nodes <- list()
  rep.fuSet.edges <- list()

  for(j in 1:result$numberOfForegrounds){
    if(!is.null(names(result$foregrounds))){
      cat(paste('INFO: processing foreground',j,'(',names(result$foregrounds)[j],')','\n'))
    }else{
      cat(paste('INFO: processing foreground',j,'\n'))
    }
    cat(paste('INFO: creating full-set result...'))
    # extracting all significantly enriched cliques
    relevantCliques <- result$relevantCliques[result$p.value[[j]][,1]<p.value]

    if(!is.null(module.reference)){
      fullSet.nodes <- sapply(get.vertex.attribute(module.reference)$name, function(id){which(database$proteins==id)})
      fullSet.edges <- get.edgelist(module.reference)
      fullSet.edges <- cbind(fullSet.nodes[fullSet.edges[,1]], fullSet.nodes[fullSet.edges[,2]])
      fullSet.edges <- t(apply(fullSet.edges,1,function(x){return(c(min(x),max(x)))}))
      fullSet.edges <- apply(fullSet.edges,1,paste,collapse="-")
      names(fullSet.nodes) <- NULL
      names(fullSet.edges) <- NULL
    }else{

    }
    cat(paste('done!','\n'))

    if(!is.null(module.reference)){
      rep.fuSet.nodes <- append(rep.fuSet.nodes, list(fullSet.nodes))
      rep.fuSet.edges <- append(rep.fuSet.edges, list(fullSet.edges))
    }

    cat(paste('INFO: creating sub-set results...','\n'))
    samples <- foreach(k = 1:repeats, .combine = 'append') %dopar% {
      runs <- unique(bootstrap.samples[[k]])
      # selecting all clique-ids that are associated with the run of that particular bootstrap sample. Additionally we restrict the search to relevantCliques only.
      ids <- intersect(relevantCliques, unique(database$results[as.numeric(database$results[,"run"]) %in% runs,"id"]))
      # interestingly, the parallelization can not access the database object since it was defined outside of this lapply. Thus, we have to store the relevant information locally...
      uniqueCliqueId <- database$uniqueCliqueId[ids]

      subSet <- do.call(rbind, lapply(uniqueCliqueId, function(clique){
        return(t(combn(unlist(strsplit(clique, split=" ")),m=2)))
      }))

      subSet.nodes <- unique(c(subSet[,1],subSet[,2]))
      subSet.edges <- unique(apply(subSet,1,function(x){paste(c(min(as.numeric(x)),max(as.numeric(x))),collapse="-")}))

      if(!is.null(module.reference)){
        diffEdges <- sum(!(fullSet.edges %in% subSet.edges))
        diffNodes <- sum(!(fullSet.nodes %in% subSet.nodes))
        diffEdges.ad <- sum(!(subSet.edges %in% fullSet.edges))
        diffNodes.ad <- sum(!(subSet.nodes %in% fullSet.nodes))
        return(c(diffNodes,diffEdges,diffNodes.ad,diffEdges.ad, list(subSet.nodes), list(subSet.edges)))
      }else{
        return(c(1,1,0,0, list(subSet.nodes), list(subSet.edges)))
      }
    }
    subSets.nodes <- lapply(seq(1,length(samples),6)+4, function(x){return(samples[[x]])})
    subSets.edges <- lapply(seq(1,length(samples),6)+5, function(x){return(samples[[x]])})

    if(is.null(module.reference)){
      pairwise.node.distance <- unlist(lapply(1:(length(subSets.nodes)-1), function(i){
        return(sapply((i+1):length(subSets.nodes), function(j){
          return(length(setdiff(union(subSets.nodes[[i]],subSets.nodes[[j]]), intersect(subSets.nodes[[i]],subSets.nodes[[j]]))))
        }))
      }))/median(sapply(subSets.nodes, length))
      pairwise.edge.distance <- unlist(lapply(1:(length(subSets.edges)-1), function(i){
        return(sapply((i+1):length(subSets.edges), function(j){
          return(length(setdiff(union(subSets.edges[[i]],subSets.edges[[j]]), intersect(subSets.edges[[i]],subSets.edges[[j]]))))
        }))
      }))/median(sapply(subSets.edges, length))
      subsampling.results <- rbind(subsampling.results, c("median-nodeStability"=1-median(pairwise.node.distance),"CB-nodeStability"=1-quantile(pairwise.node.distance, probs=0.95)[1],"median-edgeStability"=1-median(pairwise.edge.distance),"CB-edgeStability"=1-quantile(pairwise.edge.distance, probs=0.95)[1],1,1,1,1))
    }else{
      samples <- cbind(sapply(seq(1,length(samples),6)+0, function(x){return(samples[[x]])}), sapply(seq(1,length(samples),6)+1, function(x){return(samples[[x]])}), sapply(seq(1,length(samples),6)+2, function(x){return(samples[[x]])}), sapply(seq(1,length(samples),6)+3, function(x){return(samples[[x]])}))
      subsampling.results <- rbind(subsampling.results, c("mean-nodeStability"=1-median(samples[,1])/length(rep.fuSet.nodes[[j]]),"CB-nodeStability"=1-quantile(samples[,1]/length(rep.fuSet.nodes[[j]]), probs=0.95)[1],"median-edgeStability"=1-median(samples[,2])/length(rep.fuSet.edges[[j]]),"CB-edgeStability"=1-quantile(samples[,2]/length(rep.fuSet.edges[[j]]), probs=0.95)[1],"mean-additionalNodes"=median(samples[,3]),"CB-additionalNodes"=quantile(samples[,3], probs=0.05)[1],"median-additionalEdges"=median(samples[,4]),"CB-additionalEdges"=quantile(samples[,4], probs=0.05)[1]))
    }
    rep.nodes <- append(rep.nodes, list(table(unlist(subSets.nodes))/repeats))
    rep.edges <- append(rep.edges, list(table(unlist(subSets.edges))/repeats))
    cat(paste('INFO: done!','\n'))
  }

  cat("INFO: stopping parallel environment...")
  stopCluster(cl)
  cat(paste("done!",'\n'))

  rownames(subsampling.results) <- paste("Foreground-",1:result$numberOfForegrounds, sep="")
  results <- list("subSampling"=subsampling.results, "nodes"=subSets.nodes, "edges"=subSets.edges, "nodeFrequency"=rep.nodes, "edgeFrequency"=rep.edges, "fullSetNodes"=rep.fuSet.nodes, "fullSetEdges"=rep.fuSet.edges)
  return(results)
}

moduleDiscoverer.db.plotComputedPValues <- function(result=NULL, fileName=NULL, p.value=0.05){
  if(!is.null(fileName)){
    CairoPDF(file=paste(fileName, "pdf",sep="."), width=8, height=5)
  }else{
    par(mfrow=c(length(result$p.value),2))
  }

  for(i in 1:length(result$p.value)){
    if(result$computeOnTheFly){
      o <- order(result$p.value[[i]][,1], decreasing=FALSE)
      plot(result$p.value[[i]][o,2], pch=16, cex=1, type="p", col="green", main="", xlab="rank", ylab="p-value")
      points(result$p.value[[i]][o,1], pch=16, cex=1, lwd=2, type="l", lty=2, col="blue")

      plot(result$p.value[[i]][o,2], pch=16, cex=1, type="p", col="green", main="", xlab="rank", ylab="p-value", ylim=c(0,p.value), xlim=c(0, which(result$p.value[[i]][o,1]>=p.value)[1]-1))
      points(result$p.value[[i]][o,1], pch=16, cex=1, lwd=2, type="l", lty=2, col="blue")
    }else{
    }
  }

  if(!is.null(fileName)){
    dev.off()
  }
}

moduleDiscoverer.module.annotateModule <- function(module=NULL, annotation.db=NULL, annotateWith=NULL, nodeIdentifier=NULL){
  if(is.null(module)){
    cat(paste('ERROR: module is NULL!','\n'))
    stop()
  }

  if(is.null(annotation.db)){
    cat(paste('ERROR: annotation.db is NULL!','\n'))
    stop()
  }

  cat(paste("INFO: loading annotation db...",'\n'))
  tryCatch(eval(parse(text=paste("library(",annotation.db,")",sep=""))), error = function(e){cat(paste("ERROR:",e,'\n'))})
  annotation.db <- unlist(strsplit(annotation.db, split="\\."))
  annotation.db <- paste(annotation.db[-length(annotation.db)], collapse=".")

  if(is.null(annotateWith)){
    cat(paste('ERROR: annotateWith is NULL!','\n'))
    stop()
  }

  if(is.null(nodeIdentifier)){
    cat(paste('ERROR: nodeIdentifier is NULL!','\n'))
    stop()
  }

  cat(paste("INFO: preparing annotation...",'\n'))
  nodes.annotation <- get.vertex.attribute(module)$name
  annotation <- matrix(NA, ncol=length(annotateWith)+1, nrow=length(nodes.annotation))
  rownames(annotation) <- nodes.annotation
  colnames(annotation) <- c("internal",annotateWith)

  annotateWith <- paste(annotation.db, annotateWith, sep="")
  nodeIdentifier <- paste(annotation.db, nodeIdentifier, sep="")

  cat(paste("INFO: testing internal identifier...",'\n'))
  nodes.annotation <- tryCatch(eval(parse(text=paste("nodes.annotation[nodes.annotation %in% mappedkeys(revmap(",nodeIdentifier,"))]", sep=""))), error = function(e){cat(paste("ERROR:",e,"\n")); stop()})
  tmp <- tryCatch(eval(parse(text=paste("lapply(mget(nodes.annotation, revmap(",nodeIdentifier,")), paste, collapse=',')", sep=""))), error = function(e){cat(paste("ERROR:",e,"\n")); stop()})
  annotation[names(tmp),1] = unlist(tmp)

  cat(paste("INFO: collecting annotation...",'\n'))
  for(i in 1:length(annotateWith)){
    cat(paste('\t',"INFO: ",colnames(annotation)[i+1],'\n'))
    tmp2 <- tryCatch(eval(parse(text=paste("lapply(lapply(lapply(tmp, strsplit, split=','), unlist), function(x){paste(unlist(mget(x, ",annotateWith[i],")), collapse=',')})", sep=""))), error = function(e){cat(paste("ERROR:",e,"\n")); stop()})
    annotation[names(tmp2),i+1] <- unlist(tmp2)
    rm(tmp2)
  }

  cat(paste("INFO: module annotation...",'\n'))
  for(i in 1:(length(annotateWith)+1)){
    cat(paste('\t',"INFO: ",colnames(annotation)[i],'\n'))
    module <- set.vertex.attribute(graph=module, name=colnames(annotation)[i], value=annotation[get.vertex.attribute(module)$name,i])
  }

  cat(paste("INFO: done!",'\n'))
  return(module)
}

moduleDiscoverer.module.createModule <- function(result=NULL, module.name="ModuleDiscoverer - regulatory module", cores=NULL){
  if(is.null(result)){
    cat(paste('ERROR: result is NULL!','\n'))
    stop()
  }

  if(is.null(cores)){
    cores <- result$cores
  }

  if(length(result$enrichedCliques)>0){
    cat(paste("INFO: initializing parallel environment...",'\n'))
    runCores <- cores
    if(detectCores()<=cores){
      cat(paste('\t',"WARNING: Maximal number of cores detected on the system is ",detectCores(),". Using ",detectCores()-1," cores instead!",'\n',sep=""))
      runCores <- detectCores() - 1
    }
    cl <- makeCluster(runCores)
    registerDoParallel(cl)
    cat(paste("INFO: done!",'\n'))

    cat("INFO: greating graph...")
    getEdges <- function(clique){
      return(t(combn(clique, m=2)))
    }
    ec <- result$enrichedCliques
    g <- igraph::simplify(graph.edgelist(el=foreach(i = 1:length(ec), .combine = 'rbind') %dopar% {
      return(getEdges(ec[[i]]))
    }, directed=FALSE))
    cat(paste('done!','\n'))

    cat("INFO: stopping parallel environment...")
    stopCluster(cl)
    cat(paste("done!",'\n'))


    g <- set.vertex.attribute(g, "inForeground", value=get.vertex.attribute(g)$name %in% names(result$foregrounds))
    g <- set.vertex.attribute(g, "inBackground", value=get.vertex.attribute(g)$name %in% names(result$background))
    tmp <- rep(0, length(get.vertex.attribute(g)$name))
    names(tmp) <- get.vertex.attribute(g)$name
    tmp[get.vertex.attribute(g)$name %in% names(result$background)] <- 1
    tmp[get.vertex.attribute(g)$name %in% names(result$foreground)] <- 2
    g <- set.vertex.attribute(g, "ForOrBacOrNon", value=tmp)
  }else{
    g <- module.empty()
  }

  return(g)
}

moduleDiscoverer.module.performEnrichmentTest <- function(module=NULL, nodeDiseaseAssociation=NULL, identifierName=NULL, background=NULL, minAssociationsRequired=5){
  if(is.null(module)){
    cat(paste('ERROR: module is NULL','\n'))
    stop()
  }

  if(is.null(identifierName)){
    cat(paste('ERROR: identifierName is NULL','\n'))
    stop()
  }

  if(is.null(nodeDiseaseAssociation)){
    cat(paste('ERROR: nodeDiseaseAssociation is NULL','\n'))
    stop()
  }

  if(!is.data.frame(nodeDiseaseAssociation)){
    cat(paste('ERROR: nodeDiseaseAssociation has to be of type data.frame','\n'))
    stop()
  }else{
    if(dim(nodeDiseaseAssociation)[2]<2){
      cat(paste('ERROR: nodeDiseaseAssociation needs to have at least two columns named geneId and diseaseId','\n'))
      stop()
    }
    nodeDiseaseAssociation$geneId <- factor(nodeDiseaseAssociation$geneId)
    nodeDiseaseAssociation$diseaseId <- factor(nodeDiseaseAssociation$diseaseId)
  }

  cat(paste("INFO: extracting nodes...",'\n'))
  nodes = tryCatch(eval(parse(text=paste("get.vertex.attribute(module)$",identifierName,"", sep=""))), error = function(e){cat(paste("ERROR:",e,'\n')); stop()})
  nodes = nodes[nodes %in% background]
  nodes = sapply(nodes, strsplit, split=",")
  nodes = unique(nodes[nodes!=""])
  om.nodes = nodes[nodes %in% nodeDiseaseAssociation$geneId]
  om.background = background[background %in% nodeDiseaseAssociation$geneId]

  cat(paste("INFO: preparing nodeDiseaseAssociation...",'\n'))
  nodeDiseaseAssociation <- nodeDiseaseAssociation[nodeDiseaseAssociation$geneId %in% om.background,]
  nodeDiseaseAssociation$diseaseId <- factor(nodeDiseaseAssociation$diseaseId)
  nodeDiseaseAssociation <- nodeDiseaseAssociation[as.character(nodeDiseaseAssociation$diseaseId) %in% names(which(table(as.character(nodeDiseaseAssociation$diseaseId))>=minAssociationsRequired)),]
  nodeDiseaseAssociation$diseaseId <- factor(nodeDiseaseAssociation$diseaseId)

  diseaseIDs = levels(nodeDiseaseAssociation$diseaseId)
  cat(paste("INFO: testing",length(diseaseIDs),"diseases...",'\n'))
  #pb <- txtProgressBar(min = 0, max = length(diseaseIDs), style = 3)
  diseases = lapply(1:length(diseaseIDs), function(i){
    #setTxtProgressBar(pb, i)
    id = diseaseIDs[i]
    members = unique(as.character(nodeDiseaseAssociation$geneId[nodeDiseaseAssociation$diseaseId == id]))
    MD = intersect(om.nodes,members)
    MnD = setdiff(om.nodes,members)
    nMD = setdiff(members,MD)
    if(length(MD)!=0){
      return(c(id,fisher.test(matrix(c(length(MD),length(MnD),length(nMD),length(om.background)-length(MD)-length(MnD)-length(nMD)), byrow=T, ncol=2), alternative="greater")$p.value, c(length(MD),length(MnD),length(nMD),length(om.background)-length(MD)-length(MnD)-length(nMD))))
    }else{
      return(c(id,1,length(MD),length(MnD),length(nMD),length(om.background)-length(MD)-length(MnD)-length(nMD)))
    }
  })
  #close(pb)
  cat(paste("INFO: preparing results...",'\n'))
  diseases <- as.data.frame(do.call(rbind, diseases))
  colnames(diseases) <- c("diseaseId","p.value","MD","MnD","nMD","nMnD")
  diseases$p.value <- as.numeric(as.character(diseases$p.value))
  diseases <- diseases[apply(apply(as.matrix.data.frame(diseases[,c("MD","nMD")]),2,as.numeric),1,sum)>=minAssociationsRequired,]
  diseases$diseaseId <- factor(diseases$diseaseId)
  diseases <- diseases[order(diseases$p.value, decreasing=FALSE),]
  for(i in 1:dim(nodeDiseaseAssociation)[2]){
    if(!(colnames(nodeDiseaseAssociation)[i] %in% c("diseaseId","geneId"))){
      tmp <- unique(sapply(levels(diseases$diseaseId), function(id){unique(as.character(nodeDiseaseAssociation[nodeDiseaseAssociation$diseaseId==id,i]))}))
      if(length(tmp)==length(levels(diseases$diseaseId))){
        tmp <- cbind(tmp, levels(diseases$diseaseId))
        if(all(levels(diseases$diseaseId) %in% tmp[,2])){
          rownames(tmp) <- tmp[,2]
          diseases <- cbind("diseaseId"=diseases$diseaseId, tmp[as.character(diseases$diseaseId),1], diseases[,-1])
          colnames(diseases)[2] = colnames(nodeDiseaseAssociation)[i]
        }else{
          cat(paste("WARNING:",colnames(nodeDiseaseAssociation)[i],"does not match uniquely to diseaseIds --> skipping",'\n'))
        }
      }else{
        cat(paste("WARNING:",colnames(nodeDiseaseAssociation)[i],"does not match uniquely to diseaseIds --> skipping",'\n'))
      }
    }
  }
  rownames(diseases) <- diseases[,1]
  diseases <- diseases[,-1]
  cat(paste("INFO: done!",'\n'))
  return(diseases)
}

