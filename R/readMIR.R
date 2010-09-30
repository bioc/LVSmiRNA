#############################################################
##                                                         ##
## readMIR.R                                               ## 
## Author: SC                                              ##
## First version: 11 Aug 2009                              ##
##                                                         ##
## Reads in miRNA data from Agilent. Code is modified from ##
## read.maimages in the limma package                      ##
#############################################################


removeExtension <- function (x) 
{
  x <- as.character(x)
  x <- basename(x)
  
  n <- length(x)
  if (length(grep("\\.", x)) < n) 
    return(x)
  ext <- sub("(.*)\\.(.*)$", "\\2", x)
  if (all(ext[1] == ext)) 
    return(sub("(.*)\\.(.*)$", "\\1", x))
  else return(x)
}

read.mir <- function (files = NULL, path = NULL, ext = NULL, annotation=NULL,
                      names = NULL, columns = list(E="gMedianSignal"),
                      other.columns = NULL, read.bg = TRUE,
                      wt.fun = NULL, verbose = TRUE, sep = "\t", quote = "",
                      remove.ctrl=TRUE, ...)
{

  ## Background might be useful
  if(read.bg)
    {
      names(columns)[names(columns) == "E"] <- "G"
      columns <- c(columns,Gb = "gBGMedianSignal")
      
    }
  
  
  if (is.null(files))
    {
      if (is.null(ext)) 
        stop("Must specify input files")
      else
        {
          extregex <- paste("\\.", ext, "$", sep = "")
          files <- dir(path = ifelse(is.null(path), ".", path), 
                       pattern = extregex)
          files <- sub(extregex, "", files)
        }
    }
  
  
  if (is.data.frame(files))
    {
      targets <- files
      files <- files$FileName
      if (is.null(files)) 
        stop("targets frame doesn't contain FileName column")
      if (is.null(names)) 
        names <- targets$Label
    }
  else {
    targets <- NULL
  }

  slides <- as.vector(as.character(files))
  if (!is.null(ext)) 
    slides <- paste(slides, ext, sep = ".")
  nslides <- length(slides)
  if (is.null(names)) 
    names <- removeExtension(files)

##   if(missing(columns))
##     columns <- list(G = "gMedianSignal")
##   else
    if(!is.list(columns))
      stop("columns must be a named list")
  
  cnames <- names(columns)

  if(is.null(annotation))
    annotation <-  c("Row", "Col", "Start", "Sequence", "SwissProt", "GenBank", "Primate", 
                     "GenPept", "ProbeUID", "ControlType", "ProbeName", 
                     "GeneName", "SystematicName", "Description")
  
  fullname <- slides[1]
  Extn <- unlist(strsplit(basename(fullname),"\\."))
  Extn <- Extn[length(Extn)]


  funextr <- switch(Extn,"gz"=gzfile,"bz2"=bzfile,file)

  if (!is.null(path)) 
    fullname <- file.path(path, fullname)
  required.col <- unique(c(annotation, unlist(columns), other.columns))

  text.to.search <- if (is.null(wt.fun)) 
    ""
  else deparse(wt.fun)


  skip <- readHeader(fullname, funextr, columns = columns, 
                     sep = sep)$NHeaderRecords

  
  obj <- read.cols(fullname, funextr, required.col, text.to.search, 
                   skip = skip, sep = "\t", quote = "", as.is = TRUE, 
                   fill = TRUE, flush = TRUE,...)

  nspots <- nrow(obj)
  
  Y <- matrix(NA, nspots, nslides)
  colnames(Y) <- names
  RG <- columns
  for (a in cnames) RG[[a]] <- Y
  if (!is.null(wt.fun)) 
    RG$weights <- Y
    if (is.data.frame(targets)) {
      RG$targets <- targets
    }
    else {
      RG$targets <- data.frame(FileName = files, row.names = names, 
                               stringsAsFactors = FALSE)
    }
  if (!is.null(annotation)) {
    j <- match(annotation, colnames(obj), 0)
    if (any(j > 0)) 
      RG$genes <- data.frame(obj[, j, drop = FALSE], check.names = FALSE)
  }

  RG$source <- "agilent.mir"
  
  
  if (!is.null(RG$genes$Row) && !is.null(RG$genes$Col))
    {
      nr <- length(unique(RG$genes$Row))
      nc <- length(unique(RG$genes$Col))
      if (nspots == nr * nc) 
        RG$printer <- list(ngrid.r = 1, ngrid.c = 1, 
                           nspot.r = nr, nspot.c = nc)
    }
  


  if (!is.null(other.columns)) {
    other.columns <- as.character(other.columns)
    j <- match(other.columns, colnames(obj), 0)
    if (any(j > 0)) {
      other.columns <- colnames(obj)[j]
      RG$other <- list()
      for (j in other.columns) RG$other[[j]] <- Y
    }
    else {
      other.columns <- NULL
    }
  }

  for (i in 1:nslides)
    {
      if (i > 1) {
        fullname <- slides[i]
        if (!is.null(path)) 
          fullname <- file.path(path, fullname)

        obj <- read.cols(fullname, funextr, required.col, text.to.search, 
                         skip = skip, sep = "\t", as.is = TRUE, quote = "", 
                         fill = TRUE, nrows = nspots, flush = TRUE, ...)
      }

      for (a in cnames) RG[[a]][, i] <- obj[, columns[[a]]]
      if (!is.null(wt.fun)) 
        RG$weights[, i] <- wt.fun(obj)
      if (!is.null(other.columns)) 
        for (j in other.columns) {
          RG$other[[j]][, i] <- obj[, j]
        }
      if (verbose) 
        cat(paste("Read", fullname, "\n"))
    }

  if(read.bg)
    RG <- new("RGList", RG)
  else
    RG <- new("EList", RG)
  
  RG$preprocessing <- list(Background=NULL,Normalization=NULL,
                           is.log=FALSE,Summarization=NULL)
  if(remove.ctrl)
    RG <- RG[RG$genes$ControlType == 0,]
  else
    RG <- RG
  RG
  
}



readHeader <- 
  function (file, funextr, columns, sep = "\t") 
{
  if (missing(columns) || !length(columns)) 
    stop("must specify column headings to find")
  columns <- protectMetachar(as.character(columns))
  if (!length(columns)) 
    stop("column headings must be specified")
  
  con <- funextr(file, "r")
  on.exit(close(con))
  
  out <- list()

  Found <- FALSE
  i <- 0
  repeat {
    i <- i + 1
    txt <- readLines(con, n = 1)
    if (!length(txt)) 
      stop("Specified column headings not found in file")
    Found <- TRUE
    for (a in columns) Found <- Found && length(grep(a, txt))
    if (Found) 
      break
  }
  out$NHeaderRecords <- i - 1
  out$ColumnNames <- strsplit(txt, split = sep)[[1]]
  out
}


## Watchout!! Every calls to connection increases the pointer!!!
## Use of seek to place back the pointer is unsafe in Windows!!!
read.cols <- 
  function (file, funextr, required.col = NULL, text.to.search = "", sep = "\t", 
            quote = "\"", skip = 0, fill = TRUE, blank.lines.skip = TRUE, 
            comment.char = "", allowEscapes = FALSE, ...) 
{
  
  file <- funextr(file,"rt")
  on.exit(close(file))
  
  if (is.null(required.col)) 
    return(read.table(file = file, header = TRUE, check.names = FALSE, 
                      sep = sep, quote = quote, skip = skip, fill = fill, 
                      blank.lines.skip = blank.lines.skip, comment.char = comment.char, 
                      allowEscapes = allowEscapes, ...))
  text.to.search <- as.character(text.to.search)
  allcnames <- scan(file, what = "", sep = sep, quote = quote, 
                    nlines = 1, quiet = TRUE, skip = skip, strip.white = TRUE, 
                    blank.lines.skip = blank.lines.skip, comment.char = comment.char, 
                    allowEscapes = allowEscapes)
  ncn <- length(allcnames)
  colClasses <- rep("NULL", ncn)
  if (is.numeric(required.col)) {
    colClasses[required.col] <- NA
  }
  else {
    colClasses[allcnames %in% as.character(required.col)] <- NA
  }
  if (any(text.to.search != "")) 
    for (i in 1:ncn) {
      if (length(grep(protectMetachar(allcnames[i]), text.to.search))) 
        colClasses[i] <- NA
    }
  ## Don't need to skip as pointer alredy moved to the next line to be read
  secondline <- scan(file, what = "", sep = sep, quote = quote, 
                     nlines = 1, quiet = TRUE, strip.white = TRUE, 
                     blank.lines.skip = blank.lines.skip, comment.char = comment.char, 
                     allowEscapes = allowEscapes)
  if (length(secondline) > ncn) 
    colClasses <- c(NA, colClasses)

  ## go back to beginning. This might be a problem in Windows, so we check it out!!!
  if(isSeekable(file))
    invisible(seek(file,where=0,origin="start")) 
  else
    {
      close(file)
      file <- funextr(file,"rt")
    }
  read.table(file = file, header = TRUE, col.names = allcnames, 
             check.names = FALSE, colClasses = colClasses, sep = sep, 
             quote = quote, skip = skip, fill = fill, blank.lines.skip = blank.lines.skip, 
             comment.char = comment.char, allowEscapes = allowEscapes, 
             ...)
}
