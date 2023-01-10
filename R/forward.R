library("gtools")
library("matrixStats")

#' Calculating the probability of transition from multiple nodes to given node in the dag
#'
#' @param hmm Object of class List given as output by \code{\link{initHMM}},
#' @param prev_state vector containing state variable values for the previous nodes
#' @param cur_state character denoting the state variable value for current node
#' @return The Noisy_OR probability for the transition
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
#' Transprob = noisy_or(hmm=hmmA,prev_state=c("P","N"),cur_state="P") #for transition from P & N
#'                                                                    #simultaneously to P

noisy_or=function(hmm, prev_state, cur_state)
{
  l=length(which(prev_state==hmm$States[1]))
  fin=(hmm$transProbs[hmm$States[1],hmm$States[2]])^l
  fin=(1-hmm$leak_param)*fin
  if(cur_state==hmm$States[2])
    return(fin)
  else
    return(1-fin)
}

#' Calculating the probability of occurance of particular values of covariates at a node given the value of target.
#'
#' @param state character value of state variable at a particular node.
#' @param obsv character vector of values of covariates at that node.
#' @param probs emission probability distribution of the covariates in TAN structure.
#' @param pars integer vector denoting the parents of the nodes(other than root) in the TAN structure.
#' @return probability of occurance of particular values of covariates at a node given the value of target.

calc_emis=function(state, obsv, probs, pars)
{
  probab=0

  for(i in 1:length(obsv))
  {
    if(is.na(obsv[i]))
      node=" "
    else
      node=paste("\'",obsv[i],"\'",sep="")

    if(i==1)
    {
      strr=paste('pros=as.vector(probs[[',i,']][',node,',',"\'",state,"\'",'])',sep="")
      eval(parse(text=strr))
      pros=log(pros)
      probab=probab+matrixStats::logSumExp(pros)
      next
    }
    if(is.na(obsv[pars[i-1]]))
      parnode=" "
    else
      parnode=paste("\'",obsv[pars[i-1]],"\'",sep="")


    strr=paste('pros=as.vector(probs[[',i,']][',node,',',parnode,',',"\'",state,"\'",'])',sep="")
    eval(parse(text=strr))
    pros=log(pros)
    probab=probab+matrixStats::logSumExp(pros)
  }
  return(probab)
}



#' Infer the forward probabilities for all the nodes of the dagHMM
#'
#' \code{forward} calculates the forward probabilities for all the nodes
#'
#' The forward probability for state X up to observation at node k is defined as the probability of observing the sequence of observations e_1,..,e_k given that the state at node k is X.
#' That is:\cr\code{f[X,k] := Prob( X_k = X | E_1 = e_1,.., E_k = e_k)}
#' \cr where \code{E_1...E_n = e_1...e_n} is the sequence of observed emissions and \code{X_k} is a random variable that represents the state at node \code{k}
#'
#' @param hmm hmm Object of class List given as output by \code{\link{initHMM}}
#' @param observation Dataframe containing the discritized character values of only covariates at each node. Column names of dataframe should be same as the covariate names. Missing values should be denoted by "NA".
#' @param ft_seq A vector denoting the order of nodes in which the DAG should be traversed in forward direction(from roots to leaves). Output of \code{\link{fwd_seq_gen}} function.
#' @param kn_states (Optional) A (L * 2) dataframe where L is the number of training nodes where state values are known. First column should be the node number and the second column being the corresponding known state values of the nodes
#' @return (2 * N) matrix denoting the forward probabilites at each node of the dag, where "N" is the total number of nodes in the dag
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
#' ft_sq = fwd_seq_gen(hmmA)
#' kn_st = data.frame(node=c(3),state=c("P"),stringsAsFactors = FALSE)
#'                    #state at node 3 is known to be "P"
#' ForwardProbs = forward(hmm=hmmA,observation=obsvA,ft_seq=ft_sq,kn_states=kn_st)
#' @seealso \code{\link{backward}}

forward = function (hmm, observation,ft_seq, kn_states=NULL)
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
  f = array(NA, c(nStates, nObservations))
  dimnames(f) = list(states = hmm$States, index = 1:nObservations)
  for (x in 1:length(ft_seq))
  {
    k=ft_seq[x]
    bool= k %in% kn_states[,1]
    if(bool==TRUE)
    {
      ind=match(k,kn_states[,1])
      st=kn_states[ind,2]
    }
    fromstate=which(dagmat[,k]!=0)
    len_link= length(fromstate)
    if(len_link==0)
    {
      for (state in hmm$States)
      {
        f[state,k]=log(hmm$startProbs[state])
      }
      if(bool==TRUE)
      {
        st_ind=which(st!=hmm$States)
        mapdf = data.frame(old=c(1:nStates),new=hmm$States)
        tozero=as.character(mapdf$new[match(st_ind,mapdf$old)])
        f[tozero,k]=-Inf
      }
      next
    }
    prev_array=gtools::permutations(n=nStates, r=len_link, v=hmm$States, repeats.allowed = TRUE)
    inter= intersect(fromstate, kn_states[,1])
    len_inter=length(inter)
    t_value=rep(TRUE,dim(prev_array)[1])
    if(len_inter!=0)
    {
      for(i in 1:len_inter)
      {
        ind=match(inter[i],kn_states[,1])
        ind1=match(inter[i], fromstate)
        st=kn_states[ind,2]
        t_value=which(prev_array[,ind1]==st) & t_value
      }
    }
    ind_arr=which(t_value)

    for (state in hmm$States)
    {
      logsum=c()
      for (d in 1:length(ind_arr))
      {
        i=ind_arr[d]
        prev=0
        for (j in 1:dim(prev_array)[2])
        {
          prev=prev + (f[prev_array[i,j], fromstate[j]])
        }

        temp = prev + log(noisy_or(hmm,prev_array[i,],state))
        if(temp > - Inf & temp < 0)
        {
          logsum = c(logsum,temp)
        }
      }
      obsvv=observation[k,]
      emit=calc_emis(state,obsvv, hmm$emissionProbs,par_order)
      f[state, k] = matrixStats::logSumExp(logsum) + emit
    }
    if(bool==TRUE)
    {
      st_ind=which(st!=hmm$States)
      mapdf = data.frame(old=c(1:nStates),new=hmm$States)
      tozero=as.character(mapdf$new[match(st_ind,mapdf$old)])
      f[tozero,k]=-Inf
    }
  }
  return(f)
}
