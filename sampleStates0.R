cd4s <- read.csv('ukdrd2/cd4s.csv' )
res <- read.csv('ukdrd2/resistance.csv')
dem <- read.csv('ukdrd2/demographics.csv')
#~ W0 <- read.csv( '../msm0/W0.csv')
#pids <- unique( W0$donor ) #doesnt work; only ~5k in cluster are counteed here
dates <- read.table( 'ukdrd2/ExaML_trees/ExaML_result.subUKogC_noDRM.reroot_dropOG.dates', skip=1)
sids <- dates$V1
pids <- res$patientindex[ match(sids, res$testindex) ]
sampleTimes_birthday <- dates$V2

## make st rel to 1979
date0 <- as.Date('1979-01-01')
date_bd <- as.Date('1980-11-11')
sampleTimes <- as.numeric(date_bd + sampleTimes_birthday - date0 )

write.table( file = 'sampleTimes' , sampleTimes, col.names=F, row.names=F )


## sample states
# *stage, risk , care , *age 
prisk1 <- .8
dobs <- as.numeric( dem$dob_y[ match(pids, dem$patientindex) ]  )
ages <- sampleTimes / 365 + 1980 + as.numeric( as.Date('1980-11-11') - as.Date('1980-01-01')) / 365 - dobs #approximate
# age quantiles, age rates
#~ > print(qs)
#~       25%  50%  75% 100% 
#~ 18.0 27.0 33.0 40.0 80.5 
age2quantile <- function(age){
	if (is.na(age)) return (2)#NOTE should be fine for simulations, not for inference
	if (age < 27) return(1)
	if (age < 33) return(2)
	if (age < 40) return(3)
	return(4)
}
ageQuants <- sapply(ages, age2quantile)
#~ 500, 350, 200 
cd4s <- cd4s$cd4[ match(pids, cd4s$patientindex ) ]
cd4toStage <- function(cd4){
	if (is.na(cd4)) return(3) #NOTE should be fine for simulations, not for inference
	if (cd4 > 700 ) return(1) #based on .9 quantile
	if (cd4 > 500 ) return(2)
	if (cd4 > 350 ) return(3)
	if (cd4 > 200 ) return(4)
	return(5)
}
cd4Stages <- sapply( cd4s, cd4toStage )

#demes
N_NH_COMPS <- 5
N_AGE_COMPS <- 4
N_RISK_COMPS <- 2
N_CARE_COMPS <- 3
NH_COMPS <- paste(sep='', 'stage', 1:N_NH_COMPS )
AGE_COMPS <- paste(sep='', 'age', 1:N_AGE_COMPS )
RISK_COMPS <- paste( sep='', 'riskLevel', 1:N_RISK_COMPS )
CARE_COMPS <- paste(sep='', 'care', 1:N_CARE_COMPS)
DEMES <-c()
for ( nh in NH_COMPS ){
	for (age in AGE_COMPS){
		for (care in CARE_COMPS){
			for (risk in RISK_COMPS){
				DEMES <- c( DEMES, paste(sep='.', nh ,age, care, risk ))
			}
		}
	}
}
DEMES <- c( DEMES, 'src' )
m <- length(DEMES)
#pid2sampleStates <- function(pid)

sampleStates <- t( sapply( 1:length(pids), function(i) {
	nh <- paste(sep='', 'stage', cd4Stages[i])
	age <- paste(sep='', 'age', ageQuants[i])
	care <- paste(sep='', 'care', 1)
	risk <- paste(sep='', 'riskLevel', ifelse (runif(1) < prisk1, 1, 2) )
	deme <- paste(sep='.', nh ,age, care, risk )
	ss <- setNames(rep(0, m ), DEMES)
	ss[deme] <- 1
	ss
}))
write.table( file = 'sampleStates' , sampleStates, col.names=F, row.names=F )
