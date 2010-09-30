.onLoad <- function(lib, pkg){

  .C("Lapack_Initialize",PACKAGE="LVSmiRNA")

}


## .First.lib <- function(libname, pkgname) {

  
##   library.dynam("LVSmiRNA",pkgname,libname)
  
##   .C("Lapack_Initialize",PACKAGE="LVSmiRNA")

## }



.Last.lib <- function(libpath) {

  library.dynam.unload("LVSmiRNA",libpath)

}
