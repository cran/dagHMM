library(bnlearn)
library(bnclassify)

#' Generating the inital emission probability distribution of the covariates in TAN structure.
#'
#' @param net Object of type 'bn' provided as output by \link{bnlearn}{model2network} showing the TAN structure between target variable and covariates.
#' @param observation Dataframe containing the discritized character values of only covariates at each node. Column names of dataframe should be same as the covariate names. Missing values should be denoted by "NA".
#' @param sym Character vector of possible values of target variable
#' @return Inital emission probability distribution of the covariates in TAN structure
#' @examples
#'
#' library(bnlearn)
#'
#' bnet = model2network("[A][C|A:B][D|A:C][B|A]") #A is the target variable while
#'                                                #B, C and D are covariates
#' obsvA=data.frame(list(B=c("L","H","H","L","L"),C=c("H","H","L","L","H"),D=c("L","L","L","H","H")))
#' target_value= c("P","N")
#' prob= gen_emis(net=bnet,observation=obsvA,sym=target_value)


gen_emis = function(net, observation, sym)
{
  observation[]<- lapply(observation, factor)
  lpar=as.numeric(unlist(lapply(net$nodes, function (x) length(x$parents))))
  if(!(length(which(lpar==0))==1 & length(which(lpar==1))==1 & length(which(lpar==2))==(length(net$nodes)-2)))
    stop("Given structure is not TAN. Please reload the structure in valid format\n")

  target=names(net$nodes)[which(lpar==0)]
  dag_head=names(net$nodes)[which(lpar==1)]
  o_nodes=setdiff(names(net$nodes),c(target,dag_head))
  observation[[target]]=as.factor(sample(x=sym,size=dim(observation)[1],replace=T,prob = rep(1/(length(sym)),length(sym))))
  bn=bnlearn::bn.fit(net, observation, method='bayes')
  node_order=c(dag_head,o_nodes)
  emprob=list()
  for(i in 1:length(node_order))
  {
    matr=bn[[node_order[i]]]$prob
    matr[!is.na(matr)]=1/length(dimnames(matr)[[1]])
    if(i!=1)
      matr=aperm(matr,c(1,3,2))
    emprob[[node_order[i]]]=matr
  }
  return(emprob)
}


#' Initializing dagHMM with given parameters
#'
#' @param States A (2 * 1) vector with first element being discrete state value for the cases(or positive) and second element being discrete state value for the controls(or negative) for given dagHMM
#' @param dagmat Adjacent Symmetry Matrix that describes the topology of the dag
#' @param net Object of type 'bn' provided as output by \link{bnlearn}{model2network} showing the TAN structure between target variable and covariates.
#' @param observation Dataframe containing the discritized character values of covariates at each node. If "net" is not given, dataframe should also contain the column for target variable (as the last column) so as to learn the structure. Column names of dataframe should be same as the covariate names. Missing values should be denoted by "NA".
#' @param startProbs (Optional) (2 * 1) vector containing starting probabilities for the states. Default is equally probable states
#' @param transProbs (Optional) (2 * 2) matrix containing transition probabilities for the states.
#' @param leak_param (Optional) Leak parameter used in Noisy-OR algorithm used in \code{\link{forward}} and \code{\link{noisy_or}}.Default is 0
#' @return List describing the parameters of dagHMM(pi, alpha, beta, dagmat, net)
#' @examples
#'
#' library(bnlearn)
#'
#' tmat = matrix(c(0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0),
#'                5,5, byrow= TRUE ) #for "X" (5 nodes) shaped dag
#' states = c("P","N") #"P" represent cases(or positive) and "N" represent controls(or negative)
#' bnet = model2network("[A][C|A:B][D|A:C][B|A]") #A is the target variable while
#'                                                #B, C and D are covariates.
#' obsvA=data.frame(list(B=c("L","H","H","L","L"),C=c("H","H","L","L","H"),D=c("L","L","L","H","H")))
#' hmmA = initHMM(States=states, dagmat= tmat, net=bnet, observation=obsvA)
#' obsvB=data.frame(list(B=c("L","H","H","L","L"),C=c("H","H","L","L","H"),D=c("L","L","L","H","H")),
#'                       A=c("P","N","P","P","N"))
#' hmmB = initHMM(States=states, dagmat= tmat, net=NULL, observation=obsvB)

initHMM = function (States, dagmat, net=NULL,observation,startProbs = NULL, transProbs = NULL, leak_param=0)
{
  nStates = length(States)
  observation[]<- lapply(observation, factor)
  S = rep((1/nStates),nStates)
  T = 0.5 * diag(nStates) + array(0.5/(nStates), c(nStates,nStates))
  names(S) = States
  dimnames(T) = list(from = States, to = States)
  if (!is.null(startProbs)) {
    S[] = startProbs[]
  }
  if (!is.null(transProbs)) {
    T[, ] = transProbs[, ]
  }

  if(is.null(net))
  {
    TANB=bnclassify::tan_cl(colnames(observation)[dim(observation)[2]],dataset = observation)
    observation[[colnames(observation)[dim(observation)[2]]]]=NULL
    net=bnlearn::empty.graph(TANB$.dag$nodes)
    bnlearn::arcs(net)=TANB$.dag$edges
  }

  Symbols=list()
  for(i in 1:(dim(observation)[2]))
    Symbols[[i]]=as.character(sort(unique(as.character(observation[,i]))))

  E=gen_emis(net,observation,States)

  return(list(States = States, Symbols = Symbols, startProbs = S,
              transProbs = T, emissionProbs = E, adjsym=dagmat, net=net, leak_param=leak_param))
}
