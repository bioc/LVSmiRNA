\name{boxplot}

\alias{boxplot}

\title{Box Plots}

\description{

  NOTE: This man page is for the \code{boxplot} S4 generic function
  defined in this package.  See
  \code{?graphics::\link[graphics]{boxplot}} for the default methods
  (defined in the graphics package).  Bioconductor packages can define
  specific methods for objects (typically vector-like or
  data-frame-like) not supported by the default method.

}

\usage{
boxplot(x, ...)
}

\arguments{
  \item{x}{
    A vector-like or data-frame-like R object.
  }

  \item{...}{Additional arguments, perhaps used by methods.}
}

\value{

  The default \code{boxplot} method returns the value documented at
  \code{?graphics::\link[graphics]{boxplot}}. Specific methods defined
  on this generic must return similar values.

}

\seealso{
  \code{\link[graphics]{boxplot}} in the \pkg{graphics} package for the
  default \code{boxplot} methods.

  \code{\link[methods]{showMethods}} for displaying a summary of the
  methods defined for a given generic function.

  \code{\link[methods]{selectMethod}} for getting the definition of
  a specific method.

}

\examples{
boxplot
showMethods("boxplot")
selectMethod("boxplot", "ANY")  # the default method
}

\keyword{methods}
