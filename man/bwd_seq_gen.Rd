% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bwd_seq_gen.R
\name{bwd_seq_gen}
\alias{bwd_seq_gen}
\title{Calculate the order in which nodes in the dag should be traversed during the backward pass(leaves to roots)}
\usage{
bwd_seq_gen(hmm, nlevel = 100)
}
\arguments{
\item{hmm}{hmm Object of class List given as output by \code{\link{initHMM}}}

\item{nlevel}{No. of levels in the dag, if known. Default is 100}
}
\value{
Vector of length "D", where "D" is the number of nodes in the dag
}
\description{
dag is a complex graphical model where we can have multiple parents and multiple children for a node. Hence the order in which the dag should be tranversed becomes significant. Backward algorithm is a dynamic programming problem where to calculate the values at a node,
we need the values of the children nodes beforehand, which need to be traversed before this node. This algorithm outputs a possible(not unique) order of the traversal of nodes ensuring that the childrens are traversed first before the parents
}
\examples{
library(bnlearn)

tmat = matrix(c(0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0),
               5,5, byrow= TRUE ) #for "X" (5 nodes) shaped dag
states = c("P","N") #"P" represent cases(or positive) and "N" represent controls(or negative)
bnet = model2network("[A][C|A:B][D|A:C][B|A]") #A is the target variable while
                                               #B, C and D are covariates
obsvA=data.frame(list(B=c("L","H","H","L","L"),C=c("H","H","L","L","H"),D=c("L","L","L","H","H")))
hmmA = initHMM(States=states, dagmat= tmat, net=bnet, observation=obsvA)
bt_sq = bwd_seq_gen(hmmA)
}
\seealso{
\code{\link{backward}}
}
