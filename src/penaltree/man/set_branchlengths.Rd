% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo.R
\name{set_branchlengths}
\alias{set_branchlengths}
\title{Set branchlengths of tree to satisfy given nodeheights and tipheights}
\usage{
set_branchlengths(tree, nodeheights, tipheights)
}
\arguments{
\item{tree}{tree of class phylo}

\item{nodeheights}{heights of internal node}

\item{tipheights}{Height is measured as time before the most recent tip. That is, the
tip that is farthest from the root has a height of zero, the next
farthest tip has a positive height, and so on.}
}
\description{
Set branchlengths of tree to satisfy given nodeheights and tipheights
}

