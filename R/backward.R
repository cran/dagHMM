library("gtools")
library("matrixStats")

#' Infer the backward probabilities for all the nodes of the dagHMM
#'
#' \code{backward} calculates the backward probabilities for all the nodes
#'
#' The backward probability for state X and observation at node k is defined as the probability of observing the sequence of observations e_k+1, ... ,e_n under the condition that the state at node k is X.
#' That is:\cr\code{b[X,k] := Prob(E_k+1 = e_k+1, ... , E_n = e_n | X_k = X)}
#' \cr where \code{E_1...E_n = e_1...e_n} is the sequence of observed emissions and \code{X_k} is a random variable that represents the state at node \code{k}
#'
#' @param hmm hmm Object of class List given as output by \code{\link{initHMM}}
#' @param observation Dataframe containing the discritized character values of only covariates at each node. Column names of dataframe should be same as the covariate names. Missing values should be denoted by "NA".
#' @param bt_seq A vector denoting the order of nodes in which the DAG should be traversed in backward direction(from leaves to roots). Output of \code{\link{bwd_seq_gen}} function.
#' @param kn_states (Optional) A (L * 2) dataframe where L is the number of training nodes where state values are known. First column should be the node number and the second column being the corresponding known state values of the nodes
#' @return (2 * N) matrix denoting the backward probabilites at each node of the dag, where "N" is the total number of nodes in the dag
#' @examples
#'
#' library(bnlearn)
#'
#' tmat = matrix(c(0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0),
#'                5,5, byrow= TRUE ) #for "X" (5 nodes) shaped dag
#' states = c("P","N") #"P" represent cases(or positive) and "N" represent controls(or negative)
#' bnet = model2network("[A][C|A:B][D|A:C][B|A]") #A is the target variable while
#'                                                #B, C and D are covariates
#' obsvA=data.frame(list(B=c("L","H","H","L","L"),C=c("H","H","L","L","H"),D=c("L","L","L","H","H")))
#' hmmA = initHMM(States=states, dagmat= tmat, net=bnet, observation=obsvA)
#' bt_sq = bwd_seq_gen(hmmA)
#' kn_st = data.frame(node=c(3),state=c("P"),stringsAsFactors = FALSE)
#'                    #state at node 3 is known to be "P"
#' BackwardProbs = backward(hmm=hmmA,observation=obsvA,bt_seq=bt_sq,kn_states=kn_st)
#' @seealso \code{\link{forward}}

backward= function (hmm, observation,bt_seq, kn_states=NULL)
{
  if(is.null(kn_states))
    kn_states=data.frame(node=integer(),state=character(),stringsAsFactors = F)
  dagmat=hmm$adjsym
  net=hmm$net
  if(is.data.frame(observation))
  {
    observation = data.frame(lapply(observation, as.character), stringsAsFactors=FALSE)
    observation=as.matrix(observation)
    observation[which(observation=="")]=NA
  }
  target=names(net$nodes)[which(as.numeric(unlist(lapply(net$nodes, function (x) length(x$parents))))==0)]
  dag_head=names(net$nodes)[which(as.numeric(unlist(lapply(net$nodes, function (x) length(x$parents))))==1)]
  o_nodes=setdiff(names(net$nodes),c(target,dag_head))
  par_order=c()
  for(i in 1:length(o_nodes))
  {
    par2=setdiff(net$nodes[[o_nodes[i]]]$parents,c(target))
    par_order=c(par_order,which(c(dag_head,o_nodes)==par2))
  }

  hmm$transProbs[is.na(hmm$transProbs)] = 0
  nemission = length(hmm$emissionProbs)
  for(m in 1:nemission)
  {
    hmm$emissionProbs[[m]][is.na(hmm$emissionProbs[[m]])] = 0
  }
  nObservations =dim(observation)[1]
  nStates = length(hmm$States)
  b = array(NA,c(nStates, nObservations))
  dimnames(b) = list(states = hmm$States, index = 1:nObservations)
  for (x in 1:length(bt_seq))
  {
    k=bt_seq[x]
    bool= k %in% kn_states[,1]
    if(bool==TRUE)
    {
      ind=match(k,kn_states[,1])
      st=kn_states[ind,2]
    }
    nxtstate=which(dagmat[k,]!=0)
    len_link=length(nxtstate)
    if(len_link==0)
    {
      for (state in hmm$States)
      {
        b[state,k]=0
      }
      if(bool==TRUE)
      {
        st_ind=which(st!=hmm$States)
        mapdf = data.frame(old=c(1:nStates),new=hmm$States)
        tozero=as.character(mapdf$new[match(st_ind,mapdf$old)])
        b[tozero,k]=-Inf
      }
      next
    }
    next_array=gtools::permutations(n=nStates, r=len_link, v=hmm$States, repeats.allowed = TRUE )
    inter= intersect(nxtstate, kn_states[,1])
    len_inter=length(inter)
    t_value=rep(TRUE,dim(next_array)[1])
    if(len_inter!=0)
    {
      for(i in 1:len_inter)
      {
        ind=match(inter[i],kn_states[,1])
        ind1=match(inter[i], nxtstate)
        st=kn_states[ind,2]
        t_value=which(next_array[,ind1]==st) & t_value
      }
    }
    ind_arr=which(t_value)
    for (state in hmm$States)
    {
      logsum=c()
      for (d in 1:length(ind_arr))
      {
        i=ind_arr[d]
        temp=0
        for (j in 1:dim(next_array)[2])
        {
          obsvv=observation[nxtstate[j],]
          emit= calc_emis(next_array[i,j],obsvv,hmm$emissionProbs,par_order)
          temp = temp + (b[next_array[i,j], nxtstate[j]] + log(hmm$transProbs[state, next_array[i,j]]) + emit)
        }
        if(temp > - Inf & temp< 0)
        {
          logsum = c(logsum,temp)
        }
      }
      b[state, k] = matrixStats::logSumExp(logsum)
    }
    if(bool==TRUE)
    {
      st_ind=which(st!=hmm$States)
      mapdf = data.frame(old=c(1:nStates),new=hmm$States)
      tozero=as.character(mapdf$new[match(st_ind,mapdf$old)])
      b[tozero,k]=-Inf
    }
  }
  return(b)
}
