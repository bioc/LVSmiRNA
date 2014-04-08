.onLoad <- function(lib, pkg){

  .C("Lapack_Initialize",PACKAGE="LVSmiRNA")

}

.onUnload <- function(libpath) {
    
    library.dynam.unload("LVSmiRNA",libpath)
    
}
