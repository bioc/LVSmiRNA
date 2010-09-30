basicClusterInit <- function(clusterNumberNodes = 2,nameCluster="cl",
                              typeCluster = c("MPI","SOCK")) {
  ## Initialize a cluster
  ## Note that Rmpi already checks if another cluster is up and running,
  ## so we do not need to do that.

  require(snow)

  if(!(typeCluster %in% c("MPI","SOCK"))) stop("typeCluster needs to be MPI")
  ## if(!(typeCluster %in% c("MPI", "PVM"))) stop("typeCluster needs to be PVM or MPI")
  
  if(typeCluster == "MPI") {
    require(Rmpi)
  }
  ## if(typeCluster == "PVM") {
  ##   require(rpvm)
  ## }

  if(length(find(nameCluster)))
        stop("\nThere is another object called ", nameCluster,".\n",
             "It could also mean that there is\n",
             "already an object with this name; either remove the object\n",
             "or use another name for the cluster.\n")
  
  assign(nameCluster,
         makeCluster(clusterNumberNodes,
                     type = typeCluster),
         env = .GlobalEnv)    
  
  ## Load LVSmiRNA on every cluster
  clusterEvalQ(eval(parse(text = nameCluster)), library(LVSmiRNA))
}
