setwd('/work/slevu/simulations/')
require(phydynR)
source('model0.R')

MH <- 20 #10 # look up to 10 years in past for infector probs 
PID <- Sys.getpid() 

# counterfactuals sim'ed separately, eg: 
#~ nh_wtransm <- c( 
#~ 	nh1 = 1
#~ 	,nh2 = 1
#~ 	,nh3 = 1
#~ 	,nh4 = 1
#~ 	,nh5 = 1
#~ )

{

	o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
	tfgy <- .tfgy( o )

#############################
#~ incidence and prevalence 
yfin <- tfgy[[4]][[length(times_day)]]
ffin <- tfgy[[2]][[length(times_day)]]
newinf <- sum(ffin[1:120, 1:120] ) * 365 
plwhiv <- sum( yfin[-length(yfin)] )

#~  sample time and states
sampleTimes <- scan( file = 'sampleTimes' )
ss  <- matrix( scan( file = 'sampleStates' ) , byrow=TRUE, ncol = m)
colnames(ss) <- DEMES
# regularise
ss <- ss + 1e-4
ss = sampleStates <- ss / rowSums(ss)

## sim tree

print('sim tree')
print(date())
st.tree <- system.time( {
	daytree <- sim.co.tree.fgy(tfgy,  sampleTimes, sampleStates)
})
print(date())
}




# rescale tree 
tree <- daytree
sampleTimes <- days2years( tree$sampleTimes )
tree$edge.length <- tree$edge.length / 365
bdt <- DatedTree(  tree, sampleTimes, tree$sampleStates, tol = Inf)

n<- bdt$n

treeSampleStates <- tree$sampleStates
cd4s <- setNames( sapply( 1:nrow(treeSampleStates), function(k){
	deme <- DEMES[ which.max(treeSampleStates[k,] )  ]
	stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
	stage <- as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
	if (stage==1) return(1e3)
	if (stage==2) return(750)
	if (stage==3) return(400)
	if (stage==4) return(300)
	if (stage==5) return(100)
}), tree$tip.label)
ehis <-  setNames( sapply( 1:nrow( treeSampleStates), function(k){
	deme <- DEMES[ which.max(treeSampleStates[k,] )  ]
	stage <- strsplit( deme, '.' , fixed=T)[[1]][1]
	stage <- as.numeric( tail( strsplit(stage, '')[[1]], 1 ) )
	ifelse( stage==1, TRUE, FALSE) 
}), tree$tip.label)
sampleDemes <- setNames( sapply( 1:n, function(u) DEMES[which.max( treeSampleStates[u,])] ), tree$tip.label )


	W <- phylo.source.attribution.hiv( bdt
	  , bdt$sampleTimes # must use years
	  , cd4s = cd4s[bdt$tip.label] # named numeric vector, cd4 at time of sampling 
	  , ehi = ehis[bdt$tip.label] # named logical vector, may be NA, TRUE if patient sampled with early HIV infection (6 mos )
	  , numberPeopleLivingWithHIV  = plwhiv# scalar
	  , numberNewInfectionsPerYear = newinf # scalar 
	  , maxHeight = MH 
	  , res = 1e3
	  , treeErrorTol = Inf
	)  
	
	save( daytree, bdt, W, cd4s, sampleDemes, plwhiv, newinf, MH
	 , file = paste('data/simulations2/model0-simulateBaseline0/', PID, '.RData', sep=''))

