\name{@PKGNAME@}
\docType{data}
\alias{@PKGNAME@}
\title{Annotation package that collates several annotation resources.}
\description{This data object was automatically created by @AUTHOR@. It
  represents a collection of annotation packages that can be used as a
  single object named @PKGNAME@.  This object can be used with the
  standard four accessor method for all AnnotationDbi objects.  Namely:
  \code{cols}, \code{keytype}, \code{keys} and \code{select}.  Users are
  encouraged to read the vignette from the \code{OrganismDbi} package for
  more details.}
\usage{library(@PKGNAME@)}
\keyword{datasets}
\examples{
  @PKGNAME@
  cls <- columns(@PKGNAME@)
  cls
  cls <- cls[c(1,19,45)]
  kts <- keytypes(@PKGNAME@)
  kt <- kts[2]
  kts
  ks <- head(keys(@PKGNAME@, keytype=kts[2]))
  ks
  res <- select(@PKGNAME@, keys=ks, columns=cls, keytype=kt)
  head(res)
}
