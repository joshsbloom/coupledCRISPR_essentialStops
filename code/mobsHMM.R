# modified from the HMM R package
initHMM.mobs=function (States, Symbols, startProbs = NULL, transProbs = NULL, 
    emissionProbs = NULL) 
{
    nStates = length(States)
    nSymbols = length(Symbols)
    S = rep(1/nStates, nStates)
    T = 0.5 * diag(nStates) + array(0.5/(nStates), c(nStates,    nStates))
    E = array(1/(nSymbols), c(nStates, nSymbols))
    names(S) = States
    dimnames(T) = list(from = States, to = States)
    dimnames(E) = list(states = States, symbols = Symbols)
    if (!is.null(startProbs)) {
        S[] = startProbs[]
    }
    if (!is.null(transProbs)) {
        T[, ] = transProbs[, ]
    }
    if (!is.null(emissionProbs)) {
        E[, ] = emissionProbs[, ]
    }
    return(list(States = States, Symbols = Symbols, startProbs = S,     transProbs = T, emissionProbs = E))
}
 

forward.mobs=function (hmm, observation) 
{
    hmm$transProbs[is.na(hmm$transProbs)] = 0
    hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
    nObservations = length(observation)
    nStates = length(hmm$States)
    f = array(NA, c(nStates, nObservations))
    dimnames(f) = list(states = hmm$States, index = 1:nObservations)
    for (state in hmm$States) {
        f[state, 1] = sum(log(hmm$startProbs[state] * hmm$emissionProbs[state, observation[[1]]]))
    }
    for (k in 2:nObservations) {
        for (state in hmm$States) {
            logsum = -Inf
            for (previousState in hmm$States) {
                temp = f[previousState, k - 1] + log(hmm$transProbs[previousState, state])
                if (temp > -Inf) {
                  logsum = temp + log(1 + exp(logsum - temp))
                }
            }
            f[state, k] = sum(log(hmm$emissionProbs[state, observation[[k]]])) + logsum
        }
    }
    return(f)
}
backward.mobs=function (hmm, observation) 
{
    hmm$transProbs[is.na(hmm$transProbs)] = 0
    hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
    nObservations = length(observation)
    nStates = length(hmm$States)
    b = array(NA, c(nStates, nObservations))
    dimnames(b) = list(states = hmm$States, index = 1:nObservations)
    for (state in hmm$States) {
        b[state, nObservations] = log(1)
    }
    for (k in (nObservations - 1):1) {
        for (state in hmm$States) {
            logsum = -Inf
            for (nextState in hmm$States) {
                temp = b[nextState, k + 1] + sum(log(hmm$transProbs[state, nextState] * hmm$emissionProbs[nextState, observation[[k + 1]]]))
                if (temp > -Inf) {
                  logsum = temp + log(1 + exp(logsum - temp))
                }
            }
            b[state, k] = logsum
        }
    }
    return(b)
}


posterior.mobs=function (hmm, observation) 
{
    hmm$transProbs[is.na(hmm$transProbs)] = 0
    hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
    f = forward.mobs(hmm, observation)
    b = backward.mobs(hmm, observation)
    probObservations = f[1, length(observation)]
    for (i in 2:length(hmm$States)) {
        j = f[i, length(observation)]
        if (j > -Inf) {
            probObservations = j + log(1 + exp(probObservations - j))
        }
    }
    posteriorProb = exp((f + b) - probObservations)
    return(posteriorProb)
}

viterbi.mobs=function (hmm, observation) 
{
    hmm$transProbs[is.na(hmm$transProbs)] = 0
    hmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
    nObservations = length(observation)
    nStates = length(hmm$States)
    v = array(NA, c(nStates, nObservations))
    dimnames(v) = list(states = hmm$States, index = 1:nObservations)
    for (state in hmm$States) {
        v[state, 1] = sum(log(hmm$startProbs[state] * hmm$emissionProbs[state, observation[[1]]]))
    }
    for (k in 2:nObservations) {
        for (state in hmm$States) {
            maxi = NULL
            for (previousState in hmm$States) {
                temp = v[previousState, k - 1] + log(hmm$transProbs[previousState, state])
                maxi = max(maxi, temp)
            }
            v[state, k] = sum(log(hmm$emissionProbs[state, observation[[k]]])) +  maxi
        }
    }
    viterbiPath = rep(NA, nObservations)
    for (state in hmm$States) {
        if (max(v[, nObservations]) == v[state, nObservations]) {
            viterbiPath[nObservations] = state
            break
        }
    }
    for (k in (nObservations - 1):1) {
        for (state in hmm$States) {
            if (max(v[, k] + log(hmm$transProbs[, viterbiPath[k + 
                1]])) == v[state, k] + log(hmm$transProbs[state,  viterbiPath[k + 1]])) {
                viterbiPath[k] = state
                break
            }
        }
    }
    return(viterbiPath)
}


baumWelchRecursion.mobs = function(hmm, observation)
{
	TransitionMatrix    = hmm$transProbs
	TransitionMatrix[,] = 0
	EmissionMatrix      = hmm$emissionProbs
	EmissionMatrix[,]   = 0
	f = forward.mobs(hmm,  observation)
	b = backward.mobs(hmm, observation)
	probObservations = f[1,length(observation)]
	for(i in 2:length(hmm$States))
	{
		j = f[i,length(observation)]
		if(j > - Inf)
		{
			probObservations = j + log(1+exp(probObservations-j))
		}
	}
	for(x in hmm$States)
	{
		for(y in hmm$States)
		{
			temp = f[x,1] + log(hmm$transProbs[x,y]) +
					sum(log(hmm$emissionProbs[y,observation[[1+1]]])) + b[y,1+1]
			for(i in 2:(length(observation)-1))
			{
				j = f[x,i] + log(hmm$transProbs[x,y]) +
						sum(log(hmm$emissionProbs[y,observation[[i+1]]])) + b[y,i+1]
				if(j > - Inf)
				{
					temp = j + log(1+exp(temp-j))
				}
			}
			temp = exp(temp - probObservations)
			TransitionMatrix[x,y] = temp
		}
	}
	for(x in hmm$States)
	{
		for(s in hmm$Symbols)
		{
			temp = -Inf
			for(i in 1:length(observation))
			{
				if(s == observation[i])
				{
					j = f[x,i] + b[x,i]
					if(j > - Inf)
					{
						temp = j + log(1+exp(temp-j))
					}
				}
			}
			temp = exp(temp - probObservations)
			EmissionMatrix[x,s] = temp
		}
	}
	return(list(TransitionMatrix=TransitionMatrix,EmissionMatrix=EmissionMatrix))
}


baumWelch.mobs=function (hmm, observation, maxIterations = 100, delta = 1e-09, 
    pseudoCount = 0) 
{
    tempHmm = hmm
    tempHmm$transProbs[is.na(hmm$transProbs)] = 0
    tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = 0
    diff = c()
    for (i in 1:maxIterations) {
        bw = baumWelchRecursion.mobs(tempHmm, observation)
        T = bw$TransitionMatrix
        E = bw$EmissionMatrix
        T[!is.na(hmm$transProbs)] = T[!is.na(hmm$transProbs)] + 
            pseudoCount
        E[!is.na(hmm$emissionProbs)] = E[!is.na(hmm$emissionProbs)] + 
            pseudoCount
        T = (T/apply(T, 1, sum))
        E = (E/apply(E, 1, sum))
        d = sqrt(sum((tempHmm$transProbs - T)^2)) + sqrt(sum((tempHmm$emissionProbs - 
            E)^2))
        diff = c(diff, d)
        tempHmm$transProbs = T
        tempHmm$emissionProbs = E
        if (d < delta) {
            break
        }
    }
    tempHmm$transProbs[is.na(hmm$transProbs)] = NA
    tempHmm$emissionProbs[is.na(hmm$emissionProbs)] = NA
    return(list(hmm = tempHmm, difference = diff))
}












