library("future")
library("matrixStats")
library("PRROC")
library("gtools")

plan(multicore)
#' Implementation of the Baum Welch Algorithm as a special case of EM algorithm
#'
#' \code{\link{baumWelch}} recursively calls this function to give a final estimate of parameters for dag HMM
#' Uses Parallel Processing to speed up calculations for large data. Should not be used directly.
#'
#' @param hmm hmm Object of class List given as output by \code{\link{initHMM}}
#' @param observation Dataframe containing the discritized character values of only covariates at each node. Column names of dataframe should be same as the covariate names. Missing values should be denoted by "NA".
#' @param t_seq list of forward and backward DAG traversal sequence (in that order) as given output by \code{\link{fwd_seq_gen}} and \code{\link{bwd_seq_gen}} function
#' @param kn_states (Optional) A (L * 2) dataframe where L is the number of training nodes where state values are known. First column should be the node number and the second column being the corresponding known state values of the nodes
#' @param kn_verify (Optional) A (L * 2) dataframe where L is the number of validation nodes where state values are known. First column should be the node number and the second column being the corresponding known state values of the nodes
#' @return List containing estimated Transition and Emission probability matrices
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
#' kn_st = data.frame(node=c(2),state=c("P"),stringsAsFactors = FALSE)
#'                    #state at node 2 is known to be "P"
#' kn_vr = data.frame(node=c(3,4,5),state=c("P","N","P"),stringsAsFactors = FALSE)
#'                    #state at node 3,4,5 are "P","N","P" respectively
#' t_seq=list(fwd_seq_gen(hmmA),bwd_seq_gen(hmmA))
#' newparam= baumWelchRecursion(hmm=hmmA,observation=obsvA,t_seq=t_seq,
#'                              kn_states=kn_st, kn_verify=kn_vr)
#' @seealso \code{\link{baumWelch}}


baumWelchRecursion = function(hmm, observation, t_seq, kn_states=NULL, kn_verify=NULL)
{
  dagmat=hmm$adjsym
  net=hmm$net
  target=names(net$nodes)[which(as.numeric(unlist(lapply(net$nodes, function (x) length(x$parents))))==0)]
  dag_head=names(net$nodes)[which(as.numeric(unlist(lapply(net$nodes, function (x) length(x$parents))))==1)]
  o_nodes=setdiff(names(net$nodes),c(target,dag_head))
  par_order=c()
  for(i in 1:length(o_nodes))
  {
    par2=setdiff(net$nodes[[o_nodes[i]]]$parents,c(target))
    par_order=c(par_order,which(c(dag_head,o_nodes)==par2))
  }

  t_head=which(colnames(observation)==dag_head)
  o_head=which(colnames(observation) %in% o_nodes)
  observation=observation[,c(t_head,o_head)]

  observation = data.frame(lapply(observation, as.character), stringsAsFactors=FALSE)
  observation=as.matrix(observation)
  observation[which(observation=="")]=NA

  syms=hmm$Symbols
  syms=syms[c(t_head,o_head)]
  TransitionMatrix    = hmm$transProbs
  TransitionMatrix[,] = 0
  EmissionMatrix      = hmm$emissionProbs
  nemission = length(hmm$emissionProbs)


  fwd = future::future(forward(hmm,observation,t_seq[[1]],kn_states))
  bwd = future::future(backward(hmm, observation,t_seq[[2]],kn_states))
  fb_start=Sys.time()
  message("Forward_backward loop started in parallel_processes:")
  trv=FALSE
  f_count=1
  b_count=1
  while(trv==FALSE)
  {
    trv= future::resolved(fwd) & future::resolved(bwd)
    if(future::resolved(fwd)==TRUE & f_count==1)
    {
      f_time=Sys.time()
      message("Forward_loop finished")
      message('Time required: ',as.character(round(f_time-fb_start,4)),' ', units(f_time-fb_start))
      f_count=0
    }
    if(future::resolved(bwd)==TRUE & b_count==1)
    {
      b_time=Sys.time()
      message("Backward_loop finished")
      message('Time required: ',as.character(round(b_time-fb_start),4),' ', units(b_time-fb_start))
      b_count=0
    }
  }

  f=future::value(fwd)
  b=future::value(bwd)

  nStates=length(hmm$States)
  nObservations =dim(observation)[1]
  gam=f+b
  for(x in 1:length(t_seq[[1]]))
  {
    i=t_seq[[1]][x]
    summ=matrixStats::logSumExp(gam[,i])
    if(summ != -Inf)
      gam[,i]=gam[,i]-summ
  }
  if(is.null(kn_verify)==FALSE)
  {
    kn_verify[,2][which(kn_verify[,2]==hmm$States[1])]=1
    kn_verify[,2][which(kn_verify[,2]==hmm$States[2])]=0
    kn_verify[,2]=as.integer(kn_verify[,2])
    pred_prob=exp(gam["P",kn_verify[,1]])
    act_prob=kn_verify[,2]
    fg= pred_prob[which(act_prob==1)]
    bg=pred_prob[which(act_prob==0)]
    roc_obj = PRROC::roc.curve(scores.class0 = fg,scores.class1 = bg)$auc
    pr_obj = PRROC::pr.curve(scores.class0 = fg,scores.class1 = bg)$auc.integral
    message("AUC:",round(roc_obj,2))
    message("AUPR_Integral:",round(pr_obj,2))
  }

  ps_st=gtools::permutations(n=nStates,r=nStates,v=hmm$States)
  links=which(dagmat!=0,arr.ind = TRUE)
  colnames(links)=NULL
  t_prob=array(0, c(dim(links)[1],nStates,nStates))
  dimnames(t_prob)=list(c(1:dim(links)[1]),hmm$States,hmm$States)

  for(i in 1:dim(links)[1])
  {
    for(x in hmm$States)
    {
      for(y in hmm$States)
      {
        obsvv=observation[links[i,2],]
        emit=calc_emis(y, obsvv,hmm$emissionProbs, par_order)

        t_prob[i,x,y] = f[x,links[i,1]] + log(hmm$transProbs[x,y]) + b[y,links[i,2]] + emit
      }
    }
    summ=matrixStats::logSumExp(t_prob[i,,])
    if(summ!= -Inf)
      t_prob[i,,]=t_prob[i,,]-summ
  }

  for(x in hmm$States)
  {
    sumd=matrixStats::logSumExp(gam[x,t_seq[[1]]])
    for(y in hmm$States)
    {
      summ=matrixStats::logSumExp(t_prob[,x,y])
      TransitionMatrix[x,y]=exp(summ-sumd)
    }
  }
  emissmatrix=list()

  for(m in 1:nemission)
  {
    fprob=EmissionMatrix[[m]]
    pars=par_order
    sumd=matrixStats::logSumExp(gam[,t_seq[[1]]])

    if(m==1)
    {
      for(x in hmm$States)
      {
        sumd=matrixStats::logSumExp(gam[x,t_seq[[1]]])
        for(s1 in syms[[m]])
        {
          indi=intersect(which(observation[,m]==s1),t_seq[[1]])
          summ=matrixStats::logSumExp(gam[x,indi])
          fprob[s1,x] = exp(summ-sumd)
        }
      }
    }
    else
    {
      parind=pars[m-1]
      for(x in hmm$States)
      {
        for(s2 in syms[[parind]])
        {
          indp=intersect(which(observation[,parind]==s2),t_seq[[1]])
          sumd=matrixStats::logSumExp(gam[x,indp])
          for(s1 in syms[[m]])
          {
            indi=intersect(which(observation[,m]==s1),indp)
            summ=matrixStats::logSumExp(gam[x,indi])
            fprob[s1,s2,x] = exp(summ-sumd)
          }
        }
      }
    }
    emissmatrix[[colnames(observation)[m]]]=as.array(fprob)
  }
  if(is.null(kn_verify))
    return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=emissmatrix,results= gam))
  else
    return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=emissmatrix, results= list(roc_obj,pr_obj,gam)))
}




#' Inferring the parameters of a dag Hidden Markov Model via the Baum-Welch algorithm
#'
#' For an initial Hidden Markov Model (HMM) with some assumed initial parameters and a given set of observations at all the nodes of the dag, the
#' Baum-Welch algorithm infers optimal parameters to the HMM. Since the Baum-Welch algorithm is a variant of the Expectation-Maximisation algorithm, the algorithm converges to a local solution which might not be the global optimum.
#' Note that if you give the training and validation data, the function will message out AUC and AUPR values after every iteration. Also, validation data must contain more than one instance of either of the possible states
#' @param hmm hmm Object of class List given as output by \code{\link{initHMM}}
#' @param observation Dataframe containing the discritized character values of only covariates at each node. Column names of dataframe should be same as the covariate names. Missing values should be denoted by "NA".
#' @param kn_states (Optional) A (L * 2) dataframe where L is the number of training nodes where state values are known. First column should be the node number and the second column being the corresponding known state values of the nodes
#' @param kn_verify (Optional) A (L * 2) dataframe where L is the number of validation nodes where state values are known. First column should be the node number and the second column being the corresponding known state values of the nodes
#' @param maxIterations (Optional) The maximum number of iterations in the Baum-Welch algorithm. Default is 50
#' @param delta (Optional) Additional termination condition, if the transition and emission matrices converge, before reaching the maximum number of iterations (\code{maxIterations}). The difference
#' of transition and emission parameters in consecutive iterations must be smaller than \code{delta} to terminate the algorithm. Default is 1e-5
#' @param pseudoCount (Optional) Adding this amount of pseudo counts in the estimation-step of the Baum-Welch algorithm. Default is 1e-100 (Don't keep it zero to avoid numerical errors)
#' @return List of three elements, first being the infered HMM whose representation is equivalent to the representation in \code{\link{initHMM}}, second being a list of statistics of algorithm and third being the final state probability distribution at all nodes.
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
#' kn_st = data.frame(node=c(2),state=c("P"),stringsAsFactors = FALSE)
#'                    #state at node 2 is known to be "P"
#' kn_vr = data.frame(node=c(3,4,5),state=c("P","N","P"),stringsAsFactors = FALSE)
#'                    #state at node 3,4,5 are "P","N","P" respectively
#' learntHMM= baumWelch(hmm=hmmA,observation=obsvA,kn_states=kn_st, kn_verify=kn_vr)
#' @seealso \code{\link{baumWelchRecursion}}

baumWelch = function(hmm, observation,kn_states=NULL,kn_verify=NULL,maxIterations=50, delta=1E-5,pseudoCount=1e-100)
{
  t_seq=list(fwd_seq_gen(hmm),bwd_seq_gen(hmm))
  tempHmm = hmm
  nemission = length(hmm$emissionProbs)
  tempHmm$transProbs[is.na(hmm$transProbs)] = 0
  for(m in 1:nemission)
  {
    tempHmm$emissionProbs[[m]][is.na(tempHmm$emissionProbs[[m]])] = 0
  }
  diff = c()
  iter_t=c()
  auc_iter=c()
  aupr_iter=c()
  for(i in 1:maxIterations)
  {
    message("Iteration_running: ",i)
    #message('\n')
    start_time_it=Sys.time()
    # Expectation Step (Calculate expected Transitions and Emissions)
    bw = baumWelchRecursion(tempHmm, observation, t_seq, kn_states, kn_verify)
    #print(bw)
    TM  = bw$TransitionMatrix
    EM  = bw$EmissionMatrix

    # Pseudocounts
    TM[!is.na(hmm$transProbs)]    = TM[!is.na(hmm$transProbs)]    + pseudoCount
    for(m in 1:nemission)
    {
      EM[[m]][!is.na(hmm$emissionProbs[[m]])] = EM[[m]][!is.na(hmm$emissionProbs[[m]])] + pseudoCount
    }
    # Maximization Step (Maximise Log-Likelihood for Transitions and Emissions-Probabilities)
    TM = (TM/apply(TM,1,sum))
    for(m in 1:nemission)
    {
      if(length(dim(EM[[m]]))==2)
      {
        EM[[m]]=t(t(EM[[m]])/colSums(EM[[m]]))
        next
      }
      if(length(dim(EM[[m]]))==3)
      {
        for(x in hmm$States)
        {
          EM[[m]][,,x]=t(t(EM[[m]][,,x])/colSums(EM[[m]][,,x]))
        }
        next
      }
    }
    summ=0
    for(m in 1:nemission)
    {
      di= sqrt(sum((tempHmm$emissionProbs[[m]]-EM[[m]])^2))
      summ=summ+di
    }
    d = sqrt(sum((tempHmm$transProbs-TM)^2)) + summ
    message("Delta: ",round(d,10))
    diff = c(diff, d)
    tempHmm$transProbs    = TM
    tempHmm$emissionProbs = EM
    end_time_it=Sys.time()
    iter_time=difftime(end_time_it,start_time_it, units=c("auto"))
    message('Total iteration time required: ',as.character(round(end_time_it-start_time_it,4)),' ', units(end_time_it-start_time_it))
    message("\n")
    iter_t=c(iter_t,iter_time)
    if(!is.null(kn_verify))
    {
      auc_iter=c(auc_iter,bw$results[[1]])
      aupr_iter=c(aupr_iter,bw$results[[2]])
    }
    if(d<delta)
    {
      message("Convergence Reached")
      message("\n")
      break
    }
  }
  tempHmm$transProbs[is.na(hmm$transProbs)] = NA
  for(m in 1:nemission)
  {
    tempHmm$emissionProbs[[m]][is.na(hmm$emissionProbs[[m]])] = NA
  }
  if(is.null(kn_verify))
    return(list(hmm=tempHmm,stats=list(Delta_iter=diff), finprob=exp(bw$results[[3]])))
  else
    return(list(hmm=tempHmm,stats=list(Delta_iter=diff, AUC_iter=auc_iter,AUPR_iter=aupr_iter), finprob=exp(bw$results[[3]])))
}
