#~ simple model for simulating trees and testing source attribution methods
#~ Anne's cd4-staging (5 stages) X 4 age groups  X 2 levels of risk behaviour X (Undiagnosed,Diagnosed,Treated)
#~  = 120 compartments

#~ require(rcolgem)
require(phydynR) # replaces rcolgem
require(deSolve)
require(Rcpp)
sourceCpp( 'model0.cpp' )


## parameters 
age_rates <- c(agerate1 = 1/9/365
	, agerate2 = 1/6/365
	, agerate3 = 1/7/365
	, agerate4 = 1/40.5/365
) 

stage_prog_yrs <- c( .5, 3.32, 2.7, 5.50, 5.06 ) #cori AIDS
stageprog_rates <- setNames( 1 / (stage_prog_yrs  * 365 ) 
 , c('gamma1', 'gamma2', 'gamma3', 'gamma4', 'gamma5')  )


pstarts <- c( pstartstage1 = 0 #NA
 , pstartstage2 = 0.76
 , pstartstage3 = 0.19
 , pstartstage4 = 0.05
 , pstartstage5 = 0
)

theta <- c( age_assort_factor = .5
  , pRiskLevel1 = .8
  , srcMigrationRate = 1/50/365 #1 / 25 / 365 # per lineage rate of migration to source
  , srcGrowthRate = 1 / 3 / 365 # 
  , src0 = 1e3  # initial source size 
  
  , inc_scale = 0.09401734 # based on docking (see below) #.03 
  
  , max_diag_rate = 0.66227809 # based on docking (see below)  #1/3 # time 2 diag of 3 yrs
  , diag_rate_85 = 1/10 
  , accel_diag_rate = 0.03196171 # based on docking (see below ) #1 / 7 # accel of logistic function
  , treatmentEffectiveness = .95 # slows stage progression
  
	, pstarts
	, age_rates
	, stageprog_rates
)

theta_default <- theta 


nh_wtransm <- c( 
	nh1 = 1
	,nh2 = .1
	,nh3 = .1
	,nh4 = .1
	,nh5 = .3
) 

age_wtransm <- c( 
	age1 = 1
	,age2 = 1
	, age3 = 1
	, age4 = 1
)

care_wtransm <- c( 
	care1 = 1
	, care2 = .5
	, care3 = .05
)

risk_wtransm <- c( 
	risk1 = 1
	, risk2 = 10
)


## time axes & funcs
time_res <-  52 * (2013 - 1979 )  # time steps / week #1e4 #1e3 
year0 <- 1979
year1 <- 2013
date0 <- as.Date('1979-01-01')
date1 <- as.Date('2012-12-31')
times0 <- 0
times1 <- as.numeric( date1 - date0 )
times_year <- seq(year0, 2013, length.out = time_res) #to end of 2012
times_day <- seq( 0, times1, length.out = time_res )
days2years <- function( d ){
	year0  + (year1 - year0) * d / (times1 - 0 )
}
years2days <- function(y) 
{
	(times1 - times0) * (y - year0) / (year1 - year0)
}

N_NH_COMPS <- 5
N_AGE_COMPS <- 4
N_RISK_COMPS <- 2
N_CARE_COMPS <- 3
#~ also remember source 

NH_COMPS <- paste(sep='', 'stage', 1:N_NH_COMPS )
AGE_COMPS <- paste(sep='', 'age', 1:N_AGE_COMPS )
RISK_COMPS <- paste( sep='', 'riskLevel', 1:N_RISK_COMPS )
CARE_COMPS <- paste(sep='', 'care', 1:N_CARE_COMPS)

COMPS_list <- list( NH_COMPS, AGE_COMPS, CARE_COMPS, RISK_COMPS )

#~ function(i_complist, i_compPosition, deme)
#~ {
#~ ##
#~ 	if (i_compPosition == length(COMPS_list[[i_complist]])){
#~ 		j_compPosition <- 1
#~ 		j_complist <- i_complist + 1
#~ 	} else{
#~ 		j_compPosition <- j_compPosition + 1
#~ 	}
#~ 	return( list( j_complist, j_compPosition, deme ))
#~ }

NH_COORDS <- list()
AGE_COORDS <- list()
CARE_COORDS <- list()
RISK_COORDS <- list()
DEMES <-c()
k <- 1
for ( nh in NH_COMPS ){
	for (age in AGE_COMPS){
		for (care in CARE_COMPS){
			for (risk in RISK_COMPS){
				NH_COORDS[[nh]] <- c( NH_COORDS[[nh]] , k )
				AGE_COORDS[[age]] <- c( AGE_COORDS[[age]], k )
				CARE_COORDS[[care]] <- c( CARE_COORDS[[care]], k )
				RISK_COORDS[[risk]] <- c( RISK_COORDS[[risk]], k )
				DEMES <- c( DEMES, paste(sep='.', nh ,age, care, risk ))
				k <- k + 1
			}
		}
	}
}
DEMES <- c( DEMES, 'src' )
m <- length(DEMES)

# indicators for each deme; note C-indexing
NH = rep(NA, m)
AGE = rep(NA, m)
CARE = rep(NA, m )
RISK = rep(NA, m)

k <- 1
for ( care in CARE_COMPS ){
	CARE[ CARE_COORDS[[care]] ] = k -1
	k <- k + 1
}
k <- 1
for ( x in AGE_COMPS ){
	AGE[ AGE_COORDS[[x]] ] = k -1
	k <- k + 1
}
k <- 1
for ( x in NH_COMPS ){
	NH[ NH_COORDS[[x]] ] = k -1
	k <- k + 1
}
k <- 1
for ( x in RISK_COMPS ){
	RISK[ RISK_COORDS[[x]] ] = k -1
	k <- k + 1
}



# age quantiles, age rates
#~ > print(qs)
#~       25%  50%  75% 100% 
#~ 18.0 27.0 33.0 40.0 80.5 
#~ > for ( i in 2:length(qs) )
#~ + print( qs[i] - qs[i-1] )
#~ 25% 
#~   9 
#~ 50% 
#~   6 
#~ 75% 
#~   7 
#~ 100% 
#~ 40.5 






## helpers
m <- length(DEMES)
# pr row transmission goes to col
prRecipMat <- matrix( 0. , nrow = m, ncol = m )
colnames(prRecipMat) = rownames(prRecipMat) <- DEMES
.mweight <- function( rowdeme, coldeme ){
	if (rowdeme=='src') return (0)
	if (coldeme=='src') return (0)
	rowage <-  as.numeric( regmatches( rowdeme, regexec( "\\.age([0-9])", rowdeme) )[[1]][2] )
	colage <-  as.numeric( regmatches( coldeme, regexec( "\\.age([0-9])", coldeme) )[[1]][2] )
	colpss <- as.numeric( regmatches( coldeme, regexec( "stage([0-9])", coldeme) )[[1]][2] )
	colcare <- as.numeric( regmatches( coldeme, regexec( "care([0-9])", coldeme) )[[1]][2] )
	colrisk <- as.numeric( regmatches( coldeme, regexec( "riskLevel([0-9])", coldeme) )[[1]][2] )
	wcare <- ifelse( colcare == 1, 1, 0)
	wrisk <- ifelse( colrisk == 1, theta['pRiskLevel1'], 1 - theta['pRiskLevel1'] )
#~ 	browser()
	if (colpss != 1) return(0)
	wrisk * wcare  * theta['age_assort_factor']^abs( rowage - colage )
}
for (i in 1:(m-1)) for (j in 1:(m-1)){
	prRecipMat[i,j] <- .mweight( DEMES[i], DEMES[j] )
}
prRecipMat <- prRecipMat / rowSums( prRecipMat )
prRecipMat[m,] <- 0
prRecipMat[m,m] <- 1

prStageRecipMat <- matrix( 0, nrow = m, ncol = m ); 
colnames(prStageRecipMat) = rownames(prStageRecipMat) <- DEMES
.stagemweight <- function(rowdeme, coldeme){
	if (rowdeme=='src') return (0)
	if (coldeme=='src') return (0)
	rowage <-  as.numeric( regmatches( rowdeme, regexec( "\\.age([0-9])", rowdeme) )[[1]][2] )
	colage <-  as.numeric( regmatches( coldeme, regexec( "\\.age([0-9])", coldeme) )[[1]][2] )
	rowstage <- as.numeric( regmatches( rowdeme, regexec( "stage([0-9])", coldeme) )[[1]][2] )
	colstage <- as.numeric( regmatches( coldeme, regexec( "stage([0-9])", coldeme) )[[1]][2] )
	rowcare <- as.numeric( regmatches( rowdeme, regexec( "care([0-9])", coldeme) )[[1]][2] )
	colcare <- as.numeric( regmatches( coldeme, regexec( "care([0-9])", coldeme) )[[1]][2] )
	rowrisk <- as.numeric( regmatches( rowdeme, regexec( "riskLevel([0-9])", coldeme) )[[1]][2] )
	colrisk <- as.numeric( regmatches( coldeme, regexec( "riskLevel([0-9])", coldeme) )[[1]][2] )
	if ( rowstage != 1 ) return(0)
	if (colage != rowage) return(0)
	if (colcare != rowcare ) return (0)
	if (colrisk!= rowrisk) return(0)
	return( pstarts[ colstage] )
}
for (i in 1:(m-1)) for (j in 1:(m-1)){
	prStageRecipMat[i,j] <- .stagemweight( DEMES[i], DEMES[j] )
}
prStageRecipMat <- prStageRecipMat/rowSums( prStageRecipMat )
prStageRecipMat[is.na(prStageRecipMat)] <- 0


## mig mat
# NOTE uses R indices 
STAGEPROG_RECIP <- rep(-1, m)
CARE_RECIP <- rep(-1, m )
AGE_RECIP <- rep(-1, m )
#~ RISK_RECIP not needed
for (i in 1:(m-1)){
	deme <- DEMES[i] 
	age <- as.numeric( regmatches( deme, regexec( "\\.age([0-9])", deme) )[[1]][2] )
	care <- as.numeric( regmatches( deme, regexec( "care([0-9])", deme) )[[1]][2] )
	stage <- as.numeric( regmatches( deme, regexec( "stage([0-9])", deme) )[[1]][2] )
	
	if (age < length(AGE_COMPS)){
		recip_age <- age + 1
		recip_age_deme <- sub( paste(sep='', '\\.age', age)
		  , paste(sep='', '\\.age',  recip_age)
		  , deme )
		AGE_RECIP[i] = which(DEMES==recip_age_deme)
	}
	
	if (care < length(CARE_COMPS) ){
		if (care==1 || stage > 2){ # NOTE cd4 threshold for treatment
			recip_care <- care + 1
			recip_care_deme <- sub( paste(sep='', 'care', care)
			  , paste(sep='', 'care',  recip_care)
			  , deme )
			CARE_RECIP[i] = which(DEMES==recip_care_deme)
		}
	}
	
	if (stage < length(NH_COMPS) ){
		recip_stage <- stage + 1
		recip_stage_deme <- sub( paste(sep='', 'stage', stage)
		  , paste(sep='', 'stage',  recip_stage)
		  , deme )
		STAGEPROG_RECIP[i] = which(DEMES==recip_stage_deme)
	}
}


## initial conditions
y0 <- setNames( rep(0, m ), DEMES )
y0[ CARE_COORDS$care1 ] <- 1 / length( CARE_COORDS$care1 )
y0[m] <- theta['src0'] # initial source size 



#~ idea for hacking incidence and diagnosis rates(t)
#~ phillips incidence estimate -> scale so cuminf has about right value 
#~ make diagnosis rate linear from zero; tune so that 80pc diagnosed in present
# incidence (t) 
phil_inc <- data.matrix( read.table( 'incidence.tsv') )[,1] 
phil_inc_times <- seq( 1980, 2010, length.out = length(phil_inc))
d_phil_inc_times <- phil_inc_times[2] - phil_inc_times[1]
phil_inc.t <- approxfun( phil_inc_times, phil_inc, rule = 2 )
inc.t <- function(t, theta) {
	y <- days2years(t)
	i <- min(length(phil_inc), max(1, 1 + floor( (y - phil_inc_times[1]) / d_phil_inc_times )) )
	phil_inc[i] * theta['inc_scale']
}


# diagnosis rates (t)
#~ phe_diags_total <- c( 23, 	21, 	71, 	179, 	646, 	2938, 	2648, 	2385, 	1940, 	2169, 	2571, 	2847, 	2922, 	2835, 	2824, 	2931, 	2903, 	2865, 	2915, 	3267, 	3962, 	5139, 	6395, 	7409, 	7786, 	7928, 	7498, 	7388, 	7273, 	6676, 	6362, 	6219, 	6364 ) 
diag.t <- function(t, theta){
## NOTE return val needs to be rate in units of event per day 
	y <- days2years( t) 
	mdr <- theta['max_diag_rate'] #per year
	dr_accel <- theta['accel_diag_rate'] 
	dr85 <- theta[ 'diag_rate_85' ] /365 #per day
	if (y > 1985){
		return( max( dr85,  mdr / ( 1 + exp(-(y-1985) * dr_accel) ) / 365 )   ) #per day
	}
	dr85
}




# treatment rates (t)
tr.t <- function(t){
	y <- days2years(t)
	if ( y < 1995 ) return(0 )
	 ( 1 / ( 1 + exp(-(y - 2e3)/2)) ) /365 
}
#~ ys <- 1990:2012
#~ plot( ys, 1 / ( 1 + exp(-(ys - 2e3)/5)) )






## solve model 
#~ F_matrix( double incidence
#~   , NumericVector sizes
#~   , NumericVector theta
#~   , CharacterVector demes
#~   , IntegerVector NH // length m indicators for each deme 
#~   , IntegerVector AGE
#~   , IntegerVector CARE
#~   , IntegerVector RISK
#~   , NumericVector nh_wtransm // associated weight for each category
#~   , NumericVector age_wtransm
#~   , NumericVector care_wtransm
#~   , NumericVector risk_wtransm
#~   , NumericMatrix prRecipMat // pstartstage & age mixing & prisklevel
#~ )
#~ G_matrix( NumericVector sizes
#~   , NumericVector theta
#~   , CharacterVector demes
#~   , IntegerVector NH // length m indicators for each deme 
#~   , IntegerVector AGE
#~   , IntegerVector CARE
#~   , IntegerVector RISK
#~   , IntegerVector stageprog_recip // destination for migration
#~   , IntegerVector age_recip
#~   , IntegerVector care_recip
#~   , NumericVector stageprog_rates //rates for each deme
#~   , NumericVector age_rates
#~   , NumericVector care_rates // note these depend on time
#~ )

dydt <- function(t,y, parms, ... ){
	y <- pmax(y, 0 )
	incidence <- inc.t( t, theta )
	care_rates <- c( diag.t( t, theta), tr.t( t) )
	FF <- F_matrix( incidence
	  , y
	  , as.list(theta)
	  , DEMES
	  , NH
	  , AGE
	  , CARE
	  , RISK
	  , nh_wtransm
	  , age_wtransm
	  , care_wtransm
	  , risk_wtransm
	  , prRecipMat
	)
	
	GG <- G_matrix( y
	  , as.list(theta)
	  , DEMES
	  , NH
	  , AGE
	  , CARE
	  , RISK
	  , STAGEPROG_RECIP
	  , AGE_RECIP
	  , CARE_RECIP
	  , stageprog_rates
	  , age_rates
	  , care_rates
	  , prStageRecipMat
	)
	GGns <- GG
	GGns[m, ] = GG[, m ] <- 0
	
	dy <- colSums(FF) + colSums(GGns) - rowSums(GGns)
	names(dy) <- DEMES
	## aids mort 
	dy[NH_COORDS$stage5] <-  dy[NH_COORDS$stage5] -  y[ NH_COORDS$stage5 ] * stageprog_rates[5]
	## nat mort 
	dy[AGE_COORDS$age4] <-  dy[AGE_COORDS$age4] -  y[AGE_COORDS$age4] * age_rates[4]
	## source
	dy[m] <- y[m] * theta['srcGrowthRate'] 
	
#~ if (days2years( t ) > 1995-.1 ) {
#~ print( days2years( t) )
#~ csf <- colSums( FF )
#~ csg <- colSums( GGns )
#~ rsg <- rowSums( GGns )
#~ sum( csf[ CARE_COORDS$care1 ] )
#~ sum( csg[ CARE_COORDS$care1 ] )
#~ sum( rsg[ CARE_COORDS$care1 ] )
#~ sum( GGns[ CARE_COORDS$care1, CARE_COORDS$care1 ] )
#~ 	browser()
#~ }
	
	list(dy )
}


.tfgy <- function(desolve){
# for input to tree simulator 
	.t <- desolve[,1] 
	.F <- lapply( 1:nrow(desolve), function(i){
		y <- desolve[i,-1]
		t <- desolve[i, 1]
		incidence <- inc.t( t, theta )
		FF <- F_matrix( incidence
		  , y
		  , as.list(theta)
		  , DEMES
		  , NH
		  , AGE
		  , CARE
		  , RISK
		  , nh_wtransm
		  , age_wtransm
		  , care_wtransm
		  , risk_wtransm
		  , prRecipMat
		)
		rownames(FF) = colnames(FF) <- DEMES
		FF
	})
	
	.G <- lapply( 1:nrow(desolve), function(i){
		y <- desolve[i,-1]
		t <- desolve[i, 1]
		care_rates <- c( diag.t( t, theta), tr.t( t) )
		GG <- G_matrix( y
		  , as.list(theta)
		  , DEMES
		  , NH
		  , AGE
		  , CARE
		  , RISK
		  , STAGEPROG_RECIP
		  , AGE_RECIP
		  , CARE_RECIP
		  , stageprog_rates
		  , age_rates
		  , care_rates
		  , prStageRecipMat
		)
		rownames(GG) = colnames(GG) <- DEMES
		GG
	})
	
	.Y <- lapply( 1:nrow(desolve), function(i) {
		desolve[ i, -1 ] 
	})
	
	list( .t, .F, .G, .Y )
}




## unit tests
if (F)
{
t <- 1e4 
	incidence <- inc.t( t, theta )
	care_rates <- c( diag.t( t, theta), tr.t( t) )
	
	FF <- F_matrix( incidence
	  , rep(1, m )#y0
	  , as.list(theta)
	  , DEMES
	  , NH
	  , AGE
	  , CARE
	  , RISK
	  , nh_wtransm
	  , age_wtransm
	  , care_wtransm
	  , risk_wtransm
	  , prRecipMat
	)
	rownames(FF) = colnames(FF) <- DEMES 
}



if (F)
{
	tr <- sapply( times_day, function(t) diag.t( t, theta ) )
	di <- sapply( times_day, function(t) tr.t( t ) )
	inc <- sapply( times_day, function(t) inc.t( t , theta ) )
	X11(); plot( times_year, inc )
	X11(); plot( times_day, tr )
	X11(); plot( times_day, di )
}

if (F)
{
t <- 1e4 
	care_rates <- c( diag.t( t, theta), tr.t( t) )
	
	GG <- G_matrix( rep( 1, m ) 
	  , as.list(theta)
	  , DEMES
	  , NH
	  , AGE
	  , CARE
	  , RISK
	  , STAGEPROG_RECIP
	  , AGE_RECIP
	  , CARE_RECIP
	  , stageprog_rates
	  , age_rates
	  , care_rates
	  , prStageRecipMat
	)
	
	rownames(GG) = colnames(GG) <- DEMES
	
	if (F)
	{
		# NOTE this test does not work correctly...
		for (x in rownames(GG)){
			print(x )
			print(names( GG[x,which(GG[x,]  > 0)] ) )
			cat('\n\n\n' )
		}
		
		rowSums(GG)[ rowSums( GG )==0 ]
		
		for (k in 1:m){
			x <- ifelse( STAGEPROG_RECIP[k]>-1,  DEMES[STAGEPROG_RECIP[k]], NA)
			print( c( DEMES[k],x ) )
		}
		
		x <- AGE_RECIP
		x[ AGE_RECIP < 0] <- NA
		cbind( DEMES, DEMES[x] ) 
		
		x <- CARE_RECIP
		x[ CARE_RECIP < 1 ] <- NA
		cbind( DEMES, DEMES[x] )
		
		x <- STAGEPROG_RECIP
		x[ STAGEPROG_RECIP < 1 ] <- NA
		cbind( DEMES, DEMES[x] )
	}
}

#~ care 3 do not age 
#~ care 3 do not stage prog
#~ "stage3.age4.care3.riskLevel2" - no stage prog
#~ "stage5.age2.care3.riskLevel2" - no age
#~ > rowSums(GG)[ rowSums( GG )==0 ]
#~ stage2.age4.care2.riskLevel1 stage2.age4.care2.riskLevel2 
#~                            0                            0 
#~ stage2.age4.care3.riskLevel1 stage2.age4.care3.riskLevel2 
#~                            0                            0 
#~ stage5.age4.care3.riskLevel1 stage5.age4.care3.riskLevel2 
#~                            0                            0 
#~ > 


## plots
snazzy.plot <- function( desolve, agg )
{
	oa <- ( sapply( names( agg ), function( n ) {
		rowSums( desolve[, 1 + agg[[n]] ] )  ## NOTE desolve mat has time as first col
	}) )
	X11() 
	matplot( times_year, oa, type = 'l' )
	legend( x = 'topleft', legend = names(agg)
	  , col = 1:length(agg) #SPCOLS[1:length(agg)] 
	  , pch = 1
	  )
}



## debug: 
#~ o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'euler')
if (F)
{
#~ 	dydt( 0, y0, list() )
	st.o <- system.time( {
	#~ o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'adams')
	o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'euler')
	#~ o <- ode(y=y0, times=times_day, func=dydt, parms=list()  , method = 'rk4')
	})
	
	if (F){
		X11(); snazzy.plot(  o, CARE_COORDS ) 
		snazzy.plot(  o, AGE_COORDS ) 
		snazzy.plot(  o, RISK_COORDS ) 
		snazzy.plot(  o, NH_COORDS ) 
	}
	
	if (T)
	{
		tfgy <- .tfgy( o )
		
		if (F){
			n <- 12e3 #1000 #2e2
			sampleTimes <- seq( years2days( 2004 ), years2days( 2013 ), length.out = n )
			ss0 <- tfgy[[4]][[1000]]
			ss0[m] <- 0
			ss0 <- ss0 / sum(ss0)
			sampleStates <- matrix( ss0, nrow = n, ncol = m, byrow=T)
			colnames( sampleStates ) <- DEMES
		}
		
		if (T)
		{
			sampleTimes <- scan( file = 'sampleTimes' )
			ss  <- matrix( scan( file = 'sampleStates' ) , byrow=TRUE, ncol = m)
			colnames(ss) <- DEMES
			# regularise
			ss <- ss + 1e-4
			ss <- ss / rowSums(ss)
			sampleStates <- ss
			st.tree <- system.time( {
				tree <- sim.co.tree.fgy(tfgy,  sampleTimes, sampleStates
				  #, res = 1e3 #1e2#1e3
				  , step_size_multiplier= NA)
			})
			
			#plot(tree , 'phylogram', show.tip.label=F, no.margin=T, edge.width=.1, direction='downwards')
			plot(tree , 'fan', show.tip.label=F, no.margin=T, edge.width=.1, direction='downwards')
		}
		
		
		if (F)
		{
			plot(tree , 'radial', show.tip.label=F, no.margin=T)
			plot(tree , 'phylogram', show.tip.label=F, no.margin=T)
			plot(tree , 'phylogram', show.tip.label=F, no.margin=T, edge.width=.1, direction='downwards')
			plot(tree , 'radial', show.tip.label=F, no.margin=T, edge.width=.1, direction='downwards')
		}
	}
	
}


## 'dock' model
if (F)
{
	#~ PHE: 15552 diagnosed msm in london in 2012
	propDiagnosed2012 <- 4/5
	I2012 <- 15552  / propDiagnosed2012 # assuming 80pc diagnosed
	#objfun based on both the above stats: 
	objfun <- function( lntheta0 )
	{
		theta[ names(lntheta0) ] <<- exp(lntheta0) # using globals..
		o <- ode(y=y0, times=times_day, func=dydt, parms=list(), method = 'euler')
		ifin <- sum( o[nrow(o), 2:(ncol(o)-1) ] )
		idiagnosed <- sum( o[nrow(o), 1 + c( CARE_COORDS$care2, CARE_COORDS$care3) ] )
		print(paste( ifin, idiagnosed ))
		print( theta[fit_names ] )
		((ifin - I2012) / I2012)^2 + (idiagnosed/ifin - propDiagnosed2012)^2
	}
	fit_names <- c('inc_scale', 'max_diag_rate', 'accel_diag_rate')
	theta_start <- log(theta[fit_names]) # default values
	o <- optim( theta_start, objfun , control = list(trace=6, maxit=300) )
	theta_docked <- exp(o$par)
	print((o))
}


