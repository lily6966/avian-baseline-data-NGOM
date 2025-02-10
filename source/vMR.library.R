
 	
 
#-------------------------------------------------------------------------------
# GenerateRectangularPartition
#-------------------------------------------------------------------------------
# This function partitions a specified k-dimensional region into 
# an array of equally spaced k-dimensional cells with an
# option to randomize the placement of the grid cell boundaries.  
#
# DETAILS
# ----------
#	We define a parition as a set of non-overlapping 
#	hyper-rectangles or cells within some k-dimensional 
# 	cartesian coordinate system. Note, that a partition
# 	need not be exhaustive, but it is important that 
# 	cells not overlap to run with R vectorization.  
#	
#	Each cell is defined by a 
#			(cell.min, cell.max]  
#	parameter pair, one pair for each of the k dimensions.  
#	The two data.frames cell.max and cell.min 
# 	are each (C x K) dimensional matrices that contain 
#	the max/min_{i,j} values for the i-th cell, 
#		one row per cell each cell - and the j = 1...K 
#		dimesion.   
#
#	k = number of indices or dimensions used for partition.  
#		(lat, lon, time, height, group membership, etc…				   
#
#	This routine includes a flag to randomize the location 
# 	of the partition. On input you can indicate which the 
# 	the k-dimensions should be ranodmized. 
#	Randomization is accomplished by selecting a uniform
# 	random variate across cell dimenson indicated by 
# 	width parameter. On output the 
# 	vector of realized values for the k-indices
#	is produced to uniquely define the partition. 
#
# INPUT (required)
#------------------
# partition.list 
#	max = data.frame (1 x k) defining maximum of analysis extent.
#	min = data.frame (1 x k) defining maximum of analysis extent.  		   
#	width = numeric vector (k) sprecifies width of k-dimensional
# 			partition cells. 
#	randomize = logical vector(k) TRUE or FALSE indicating whether to randomize
#		cell location in given dimension 
#
# INPUT (optional)
#------------------
# rng.seed = NULL: integer value used to control 
# 		random number generation. 
# 			 
# OUTPUT
#--------
# partition.list 
#	max = data.frame (nc x k) 
#		nc= number of cells in partition. 
#		each value defines the maximum for j-th dimension of 
#		i-th cell in paratition. j = 1 … k. i = 1….nc.
#		Same defintion for min.
#	min =  (nc x k) nc= number of cells in partition. 
#	inits = vector of realized values for the k-indices
#		used to randomize or "jitter" the partition. 
# 		Will return zeros if randomize = F. 
#
# NOTES & TODO
# -------------
#** Check inputs to make sure input dimensions match
#** ASSUMED THAT INDEXES ARE IN THE SAME COLUMN ORDER!!!!
#** Add logic so that index columns are matched by column names
#	instead of column order.
#** I am assuming that all partition indices are continuous 
#	number indices. NO FACTORS!! Can/Should we extend this to factors? 
#** names(partition.list$max) == assumed names for partition indices
#** Is there an easy way to extend this so that I can "wrap" an index? 
#	Eg day of the year for modeling interannual migration?  
#** THEO - Check out quadtree function in Matlab. 
#		is it available in R?
#** NOTE: In many applications many of the partition cells are 
# 	"empty". The strategy used here is to let the user 
# 	decide how to decide/define when a given cell is empty. 
#	The user can control this in the user-defined model training 
#	and model prediction functions. They can specify how 
# 	much data is needed within a cell before training or 
# 	prediction takes place. 
#-------------------------------------------------------------------------------
GenerateRectangularPartition <- function(
	partition.list,
	rng.seed = NULL){
#-------------------------------------------------------------------------------
if (!is.null(rng.seed)) set.seed(rng.seed)
k <- ncol(partition.list$max)
u.inits <- rep(0, k)
grids.min <- vector("list", k)
grids.max <- vector("list", k)
for (iii in 1:k){
  if (partition.list$randomize[iii]){
  	u.inits[iii] <- runif(n=1, min=0, max=1) * partition.list$width[iii]
  	} 		
  ttt <-  c(partition.list$min[1,iii], 
    		seq(from= partition.list$min[1,iii] + u.inits[iii], 
				to =  partition.list$max[1,iii],
				by =  partition.list$width[iii]),
			partition.list$max[1,iii])
	ttt <- sort(unique(ttt))
	grids.min[[iii]] <- ttt[1:(length(ttt)-1)]
	grids.max[[iii]] <- ttt[2:length(ttt)]
	rm(ttt)
	}	
min.vals <- expand.grid(grids.min)
names(min.vals) <- names(partition.list$max)
max.vals <- expand.grid(grids.max)
names(max.vals) <- names(partition.list$max)
return(list(
		max=max.vals, 
		min=min.vals, 
		inits=u.inits))
#-------------------------------------------------------------------------------
} # END FUNCTION 
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# RectangularPartitioning
#-------------------------------------------------------------------------------
# This function assigns each data point to at most a
# single cell within the given partition. 
#
# DETAILS
# ---------
#	NOTE - The map.parameters indices and 
#	partition.data indices must have exactly the same column names.  
#
# INPUT
# ------
#	partition.data = (n x k) data.frame of the data locations.
# 	map.parameters = list object with these components	
# 		max = data.frame (C x k) 
#		min = data.frame (C x k) defining max or min of extent 
# 			in k dimensions. One cell per row.
#
# OUTPUT
# -------- 
#	cell.membership = (n x 1) data.frame of cell membership 
# 		for all data pts. Cell numbers correspond to the row
# 		number in the partition.list. Any data point not assigned 
#		to a cell in the partition returns NA. 
#
# TODO
# -----------
#** Check inputs to make sure partition.list and data.index dimensions match
# 	and names match. 
#** Rewrite with apply() & or outer() to avoid the loops
# 	Consider if this route is scalable in terms of memory. 
# 	because outer produce of n locations X C cells is a potentially 
# 	very large matrix!
#-------------------------------------------------------------------------------
RectangularPartitioning <- function(
	partition.data, 
	map.parameters){ 
#-------------------------------------------------------------------------------
k <- ncol(map.parameters$min)
C <- nrow(map.parameters$min)
n <- nrow(partition.data)
cell.membership <- rep(NA, n) 
# Match indices using column names - assumes identical names
# but order does not matter, and extra columns are OK. 
iii.min.name <- rep(NA,k)
iii.max.name <- rep(NA,k)
for (iii in 1:k){
	iii.min.name[iii] <- c(1:ncol(partition.data))[
			names(partition.data) == names(map.parameters$min)[iii]   ]
	iii.max.name[iii] <- c(1:ncol(partition.data))[
			names(partition.data) == names(map.parameters$max)[iii]   ]
}
for (jjj in 1:C){
	ttt <- rep(T, n)
	for (iii in 1:k){
	  # Match Indices using column order - assumes identical order
	  # ttt <- ttt & (partition.data[,iii]  > map.parameters$min[jjj,iii] & 
	  # 			partition.data[,iii] <= map.parameters$max[jjj,iii]) 
	  # Match indices using column names	  
	  ttt <- ttt & (partition.data[,iii.min.name[iii]]  > 
	  					map.parameters$min[jjj, iii  ] &
	  				partition.data[,iii.max.name[iii]]  <= 
	  					map.parameters$max[jjj, iii ] )
	}#k
	cell.membership[ttt] <- jjj
	rm(ttt)
}
return(cell.membership)
#-------------------------------------------------------------------------------
} # END FUNCTION 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot Rectangular Partition, 2D
#-------------------------------------------------------------------------------
PlotRectangularPartition <- function(
	partition.list, x.name, y.name, ... ){
#-------------------------------------------------------------------------------
k <- ncol(partition.list$min)
C <- nrow(partition.list$min)
#  Match indices using column names - assumes identical names
x.min <- partition.list$min[, x.name == names(partition.list$min)]
x.max <- partition.list$max[, x.name == names(partition.list$max)]
y.min <- partition.list$min[, y.name == names(partition.list$min)]
y.max <- partition.list$max[, y.name == names(partition.list$max)]
# --------
plot(c(x.min,x.max),
	 c(y.min,y.max),
	type="n", ... )
for(jjj in 1:C)
	polygon(c(x.min[jjj], x.min[jjj], 
			x.max[jjj], x.max[jjj]),
		c(y.min[jjj], y.max[jjj],
			y.max[jjj],y.min[jjj]))
} # END FUNCTION
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# THE ENSEMBLE LEVEL FUNCTIONS
#-------------------------------------------------------------------------------
# Ensemble level - one goal is to hide all of the 
# details about the ensemble of partitions and all of 
# the detail about cells within a partition. 
#
# All ensemble Results need to be indexed - 
# either by the ensemble and cell number or by the 
# an observation data location.  
#
# I have broken up the computation of training, prediciton, and 
# other model diagnostics into three separate MR-like operations. 
# This is because each of these three tasks has distinct 
# input parameters and output types/products 
# For model training, the predictive model needs
# 		The training data (responses & predictors, possible some weights), and 
# 		model parameters. 
# The essential output is the model object itself, 
# which are saved as list. 
#
# For prediction, each basemodel needs 
# 		prediction data frame (predictors)  (m x p) 
# 		and the list of model objects
# and will produce a single numeric response for each 
# of the m predicted observations. The predictions 
# need to be returned as a (m x 1) vector with the 
# row order corresponding too the prediction data. 
#
# There are other model diagnostics that may require 
# addition parameters on input and may produce different 
# output types. For example, if we want to collect 
# statistics on the relative importance of the p 
# predictors used in the model, then the output 
# will be an (m x p) matrix or data frame. 
# I have not done this yet. 
# 
# TODO: 
# -------
# 	multiple randomized partitions
# 	and resampling (Friedman type justification?)
#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Ensemble of Partitions - Model Training
#-------------------------------------------------------------------------------	
# The code is a wrapper around vMR that creates and collects results 
# from baseModel training across an ensemble of partitions. 
#
# INPUTS
# -------
# Define Partition Ensemble
# 	n = number of partitions in ensemble
# 	generate.ensemble = the function to generate ensemble of partitions
# 	ensemble.parameters = parameters for generate.ensemble, used 
#		to intialize the ensemble when randomized.   
#
# Define Parititon - the Map Step
# 	partition.function = the partition function is assumed to have 
# 		two parameters partition.data & partition.parameters. 
#			partition.function(partition.data, partition.parameters)
# 		In this routine, the partition.parameters are generated 
# 		internally via the generate.ensemble function 
# 	partition.data = (n x k) data.frame: data.index object
#
# Within Partition Reduce or Aggregate step
#	baseModel.train, e.g. GAM_BaseModelTrain	
# 		This is the "reduce" or aggregation step. I think of 
# 		model training as a smoother that essentially produces
# 		a generalization of the observed data with reduced dimension. 
#	baseModel.data = (n x p) data.frame 
#		Entry X_i,j of the data frame contains the the 
#		j=1,...,p covariate values for the i-th observation
# 		where i = 1,…,n and n = sample size.
#	... = Pass through baseModel parameters
# 	baseModel.vectorize.parameters = NULL 
#		user has control over which basemodel parameters are vectorized- 
#		ie which parameters are local-call and which are global
#			list of parameters
# 		passed to BaseModelTrainModelFunction(). 
# 		Note that this offers a lot of flexibility. 
# 		using parameters to adaptively train the base model. 
# 		other parameters could be passed all the way to the 
# 		base model 
#
# OUTPUTS
# -----------------------
# list with two components. 
#	partition = list of length n. partition[[i]] must contain 
# 		all of the information to define the i-th partition 
# 		of the ensemble.   
#	model = list of length n. model[[i]] is the LIST of model objects
# 		returned for the i-th partition of the ensemble. 
# 	   
#
# TODO
# -----------------------
# Data Resampling Init Parameters
# 	Step allowing data resampling for each partition 
# 		unique-locations
# 		checkerboard 
# 		maximal coverage - equal weight
# 		quad-trees 
#
#-------------------------------------------------------------------------------
EnsembleModelTrain <- function( 
	n, 
	generate.ensemble, 
	ensemble.parameters, 
	partition.function,
	partition.data, 
	baseModel.train, 
	baseModel.data,
	...,
	baseModel.vectorize.parameters=NULL){
#-------------------------------------------------------------------------------
ens.partition <- vector("list", n)
ens.model <- vector("list", n)
vectorize.args <- NULL
if (!is.null(baseModel.vectorize.parameters))
		vectorize.args <- c(baseModel.vectorize.parameters)
for (iii in 1:n){	
	ens.partition[[iii]] <- generate.ensemble( ensemble.parameters )
	ens.model[[iii]] <-  vMR(
			map.function = partition.function, 
			map.data = partition.data, 
			map.parameters = ens.partition[[iii]],
			reduce.function= baseModel.train,
			reduce.data = baseModel.data,
			...,
			reduce.parameter.vectorize = vectorize.args)	 	
}
return(list(partition = ens.partition, 
			model = ens.model))
#-------------------------------------------------------------------------------
}# END FUNCTION
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# subsampling.function
# 
# This function is for simple subsubsampling. 
# NOTE, this function will not support 
# sampling with replacement because of the way the sample is constructed 
# when this return vector is passed back to EnsembleModelTrainSS.
# 
# INPUT: 
# 	partition.data - n x k - data frame of data locations
# 	resampling.parameters$fraction - the fraction to sample. 
# OUTPUT:
# 	index - logical vector of length nrow(partition.data)
# 			TRUE indicates selected sample
#-------------------------------------------------------------------------------
simple.sampler <- function( 
		partition.data, 
		parameter.list){
n <- nrow(partition.data)
return.vector <- rep(F, n)
# -------------
if (parameter.list$fraction > 1)  parameter.list$fraction <- 1
if (parameter.list$fraction <= 0) 
	stop("ERROR: simple.sampler() - parameter.list$fraction <= 0")
ssize <- round( parameter.list$fraction * n)
if (ssize == 0) 
	stop("ERROR: simple.sampler() - round( parameter.list$fraction * n) == 0")
# -------------
return.vector[ sample( c(1:n), size=ssize) ] <- T
return(return.vector)
} # END FUNCTION
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Ensemble Training with Partition-Level SUBSAMPLING
#
# I want to include functionality for data resampling
# across partitions to enable better estimates 
# predicted estimates and better estimates of 
# prediction variability. This routine extends the 
# last version by adding a resampling function. 
# 
# For very large data problems, it may be more 
# efficient to resample data on a cell by cell basis.
# In this case the resampling can be included into 
# the analysis by the user in the baseModel functions. 
#
# ------------------------------------------
#
# Additional OUTPUT: 
# -------------------
#  	sample = list of length n. sample[[i]] 
# 		is a logical vector with length equal to the number of
# 		rows in the training data. 
#
# TO DO: 
# --------
# change parameter to fract of training sample 
#-------------------------------------------------------------------------------
EnsembleModelTrainSS <- function( 
	n, 
	generate.ensemble, 
	ensemble.parameters, 
	partition.function,
	partition.data, 
	baseModel.train, 
	baseModel.data,
	...,
	baseModel.vectorize.parameters = NULL,
	resampling.function = NULL,
	resampling.parameters = NULL ){
#-------------------------------------------------------------------------------
ens.partition <- vector("list", n)
ens.model <- vector("list", n)
if (!is.null(resampling.function))
	ens.sample <- vector("list", n)
vectorize.args <- NULL
if (!is.null(baseModel.vectorize.parameters))
		vectorize.args <- c(baseModel.vectorize.parameters)
for (iii in 1:n){	
	ens.partition[[iii]] <- generate.ensemble( ensemble.parameters )
	ss.index <- c(1:nrow(partition.data))
	if (!is.null(resampling.function)){
		# ss.index - logical vector of length nrow(partition.data)
		ss.index <- resampling.function(partition.data,
										resampling.parameters )
		ens.sample[[iii]] <- ss.index
	}		
	# browser()
	# if partition.data has only one index, or column, 
	# then you need to manually add the col names
	# because taking the ss.index subset returns a numeric!!
	# you have to love R.
	par.data <- as.data.frame(partition.data[ss.index, ])
	names(par.data) <- names(partition.data)			
	ens.model[[iii]] <-  vMR(
			map.function = partition.function, 
			map.data = par.data, 
			map.parameters = ens.partition[[iii]],
			reduce.function= baseModel.train,
			reduce.data = baseModel.data[ss.index, ],
			...,
			reduce.parameter.vectorize = vectorize.args)	 	
	rm(ss.index,par.data)
}
if (!is.null(resampling.function))
	return.list <- list(partition = ens.partition, 
			model = ens.model,
			sample = ens.sample)
if (is.null(resampling.function))
	return.list <- list(partition = ens.partition, 
			model = ens.model)
return(return.list)
#-------------------------------------------------------------------------------
}# END FUNCTION
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# stem & predict.stem
# 
# These are wrappers for EnsembleModelTrainSS & . 
# The idea is to simplify the call by setting reasonable defaults
# and obscuring details and extra options.   
#
# NOTE- It is assumed that every dimension of the location 
# axes are randomized.  
#
#
#-------------------------------------------------------------------------------
stem <- function( 
	data, 		# train.data
	locations, 	#train.locations 					locs = train.locations,
	width,  	# = c(0.25, 0.25, 0.1),
	n,			#n = nnn.ens,
	baseModel,  #GAM_train, 
	... ,		# pass parameters to base model
	fraction = 1.0){
#-------------------------------------------------------------------------------
ttt.min <- apply(locations, 2, min)	
ttt.df <- matrix(ttt.min,1,ncol(locations))
ttt.df <- as.data.frame(ttt.df)
names(ttt.df) <- names(ttt.min)	
ttt.max <- apply(locations, 2, max)	
ttt.max.df <- matrix(ttt.max,1,ncol(locations))
ttt.max.df <- as.data.frame(ttt.max.df)
names(ttt.max.df) <- names(ttt.max)
# ---------	
init.partition <- list( 
	min = ttt.df ,
	max = ttt.max.df, 
	width = width,
	randomize = rep(T,ncol(ttt.df)))
# ---------
return.list <- EnsembleModelTrainSS( 
	baseModel.data = data,
	partition.data = locations,
	n = n,
	baseModel.train = baseModel, 
	... ,
	# Generate Ensemble of partitions
	generate.ensemble = GenerateRectangularPartition, 
	ensemble.parameters = init.partition, 
	# Partition Resampling Parameters
	resampling.function = simple.sampler,
	resampling.parameters = list(fraction = fraction),   
	# Map or Partition Step
	partition.function = RectangularPartitioning)
return(return.list)
#-------------------------------------------------------------------------------
}# END FUNCTION
#-------------------------------------------------------------------------------
predict.stem <- function(
			stem.obj, 
			new.data, 
			new.locations,
			baseModel,
			... ){
# ------------------------------------------ 
pred.stem <- EnsembleModelPredict(  
		partition.function = RectangularPartitioning,
		partition.data = new.locations, 
		ensembleModel = stem.obj,
		baseModel.predict = baseModel,  
		baseModel.data = new.data,
		...)	
return(pred.stem)
#-------------------------------------------------------------------------------
}# END FUNCTION
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Exsemble of Partitions - Model Prediction
#------------------------------------------------------------------------------- 
# Ensemble predictions are calculated as the average across all 
# parititions, drilling down at the same index location. 
# For each prediction index (row) this function 
# will return the average prediction and SD across 
# all partition predictions. 
#
# Predictions are a fundamental quantity becasue they can 
# be used to systematically characterize a model and
# with additional assumptions, make statistical inferences.
# See Hastie, Tibshirani, & Friedman.   
#
# INPUTS
# -------
# 
# Define Parititon - the Map Step
# 	partition.function 
# 	partition.data = (n x k) data.frame: data.index object
# EnsembleModel = list output from EnsembelmodelTrain
# 		defines ensemble of partitions and the models.
# Within Partition Reduce or Aggregate step
#	baseModel.predict	
#	baseModel.data = (n x p) data.frame 
#		Entry X_i,j of the data frame contains the the 
#		j=1,...,p covariate values for the i-th observation
# 		where i = 1,…,n and n = sample size.
#	... = Pass through baseModel parameters
# 	baseModel.vectorize.parameters = NULL 
#		user has control over which basemodel parameters are vectorized- 
#		ie which parameters are local-call and which are global
#			list of parameters
# 		passed to BaseModelTrainModelFunction(). 
# 		Note that this offers a lot of flexibility. 
# 		using parameters to adaptively train the base model. 
# 		other parameters could be passed all the way to the 
# 		base model 
# 	ensmeble.members = NULL, or a numeric vector of a subset of 
# 		predictions & models over which to make and sum 
# 		ensemble predictions. 
#
# OUTPUTS
# -----------------------
# Predictions are made cell by cell and then 
# repacked into the original row order of the 
# prediction design, saving three sufficient statistics:
# 	mean
#	variance
#	and sample size
# rows for which there was no prediction are returned 
# with NA's. NA's can arise for two different reasons. 
# Either that data location was not in any cell in 
# the partition or the data was located in the cell 
# but the base model did not make the prediction. 
# 	   
#
# TODO
# -----------------------
# Give user ability to collect quantiles across ensemble. 
#  
# We have to include logic that will check all of the 
# vectorized parameters to make sure that they have length = 
# to length(split.data). Then we will need to assmeble the corresponding 
# data object, like I do with ttt.mo
# the temporary list of Model Objects 
#
#-------------------------------------------------------------------------------
EnsembleModelPredict <- function(  
	partition.function,
	partition.data, 
	ensembleModel,
	baseModel.predict,  
	baseModel.data,
	...,
	baseModel.vectorize.parameters = NULL,
	ensemble.members = NULL ){
#------------------------------------------------------------------------------- 
vectorize.args <- c("model.obj")
if (!is.null(baseModel.vectorize.parameters))
	vectorize.args <- c("model.obj",baseModel.vectorize.parameters)
if (is.null(ensemble.members))
	ensemble.members <- c(1:length(ensembleModel$model))
ens.count <- rep(0, nrow(baseModel.data))
ens.mean <- rep(0, nrow(baseModel.data))
ens.var <- rep(0, nrow(baseModel.data))
# ----------------------------------------
for (jjj in ensemble.members){	
	map <- partition.function(partition.data, ensembleModel$partition[[jjj]])
    split.data.index <- split(c(1:nrow(baseModel.data)), map)
 	# ---------------------------------------------------------
	# Prediction locations may come from only
	# a subset of the cells in the partition, or none 
	# at all. Its essential that the vectorized function
	# ALL of the vectorized arguments 
	# 		data (i.e. split data) and the model object list
	# have the same number of 
	# objects in their list and that the order corsesponds.
	#
 	# The prediction data and the models intersect 
 	# with possibily different cells in the partition. 
 	# We identify the cells in common and construct 
 	# the data list and model list corresponding 
 	# to these cells. 	
 	#
 	# Go through every data.cell cell by cell
 	# if (model exists for this cell)
 	# 	add model to model.list
 	# if (model does not exist for this cell, 
 	# 		then enter NA for model.list
 	# ---------------------------------------------------------
 	data.cells  <- as.numeric(names(split.data.index))
 	if ( length(data.cells) > 0) {
 		model.cells <- as.numeric(names(ensembleModel$model[[jjj]])) 		
		model.list <- vector("list", length(data.cells))
		for (kkk in 1:length(data.cells)){
			model.list.index <- model.cells == data.cells[kkk] 
			model.list[[kkk]] <- NA
			if (sum(model.list.index)>0) {
				model.list.index <-c(1:length(model.cells))[model.list.index]
				model.list[[kkk]] <- ensembleModel$model[[jjj]][[ 
							model.list.index   ]]
			}
		} #kkk
	ens.pred <-  vMR(
			map.function = partition.function, 
			map.data = partition.data, 
			map.parameters = ensembleModel$partition[[jjj]],
			reduce.function= baseModel.predict,
			reduce.data = baseModel.data,
			model.obj = model.list,
			reduce.parameter.vectorize = vectorize.args,
			...)
	# ens.pred is a list of predicitons
	# that need to be repacked into the right order. 
	# Each list's objects are returned in same row order 
	# as the subset of the data passed on input.
	for (iii in 1:length(ens.pred)){  
		# Check if model was not fit 
		if (!is.null( ens.pred[[iii]] ) ) { 
			# Accumulate Predictions
			partition.index <- split.data.index[[iii]]		
			ens.count[partition.index] <- ens.count[partition.index] + 1
			ens.mean[partition.index] <- ens.mean[partition.index] + 
												ens.pred[[iii]] 
			ens.var[partition.index] <- ens.var[partition.index] + 
												ens.pred[[iii]] ^2
			} # if 
		} # iii 
	rm(ens.pred, model.list, partition.index)	
 	} # end length data.cells > 0 							 	
} # jjj
# ----------------------------------------------------------
# Return mean and variance
# -----------------------------------------------------------
na.ind <- (ens.count > 0)
ens.mean[na.ind] <- ens.mean[na.ind]/ens.count[na.ind]
ens.mean[!na.ind] <- NA
na.ind2 <- (ens.count >= 2)
ens.var[na.ind2] <- ens.var[na.ind2]/ (ens.count[na.ind2]-1) -
						(ens.mean[na.ind2]^2)*(ens.count[na.ind2])/
						(ens.count[na.ind2]-1)
ens.var[!na.ind2] <- NA
# -------------------------------------
return(list(mean=ens.mean, 
			var=ens.var, 
			count=ens.count))
#-------------------------------------------------------------------------------
}# END FUNCTION
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# vMR 
#-------------------------------------------------------------------------------
# MapReduce engine for Predictive Models. I have made limited modifications
# to the mapReduce() function & package by adding a little more 
# flexibility to pass parameters to the map and reduce 
# functions. 
#
# INPUTS: 
# ---------
# MAP INPUTS
# 	map.function = The map.function partitions the map data 
# 					according to the map parameters
#	map.data = (n x k) data.frame - a data index object with location variabes
#				n = sample size
# 				k = dimension of indices
#	map.parameters = partition.list object = list object with these components	
# 		max = data.frame (C x k) 
#		min = data.frame (C x k) defines max or min of extent in k dimensions. 
#					One cell per row. C = number of cells in partition
# REDUCE INPUTS
#	reduce.function = e.g. GAM_predict
#	reduce.data  = this is the the data object passed that is "reduced" 
# 					or aggregated within each partition cell.
# 					For predictive model training, this is the 
#					data design (Predictor design + Response variable)
#					The BaseModel functions determines how exactly 
# 					this data frame is used. Eg. whether all 
# 					columns/predictors are used, or some are used, etc. 
# 					 
# OPTIONAL REDUCE INPUTS
# 	... The elipses are the reduce.function parameters that 
# 		are passed through to the BaseModel fucntion.  
# 		This makes vMR very fleixble, because the user is given 
# 		control of passing parameters to the BaseModel functions. 
# 		It is possible to pass both 
# 		global and locally varying (vectorized) parameter values. 
#			(in development for local parameters, 
# 			this will require a little more logic in this routine, 
#			see the TODO list)  
#	reduce.parameter.vectorize = NULL,	
# 		This is a vector of the parameters names to be vectorized. 
# 		For example, for predictive modeling done here the 
# 		data object "D" is vectorized. By vectorized, we mean 
# 		that the "D" passed into is the data for THAT partition cell! 
# 		For doing predictions and diagnostics, the input model 
# 		object "reduce.model.obj" is also vectorized. 
#
# OUTPUT
# --------
#  vMR_result = list of C results, one for each cell in the partition. 
# 			NA is returned for cells where reduce fucntion was 
# 			not applied. 
#
# DETAILS
# ------------
# The routine works by creating a vectorized version of the 
# reduce function. One cool aspect of this is that 
# the user defines which of the parameters are also 
# vectorized. Then use ... to pass through parameters to the 
# reduce function. 
# 
# Map Function Requriements
# --------------------------------  
# 	map <- map.function(map.data, map.parameters)
# 	Takes two arguments, the data and paramters. 
# 	NA is returned for data points that are not assigned to 
# 	any of the partition cells. 
#
# Reduce Function Requirements
# --------------------------------   
# 	The reduce function must include the 
# 	arguments "D" plus any other arguments named in
# 	reduce.parameter.vectorize. 
# 	Returns NA when not applied - e.g. when minimum cell 
# 	sample size is not met. 
#
# SOURCE PROVENANCE & EVOLUTION 
# ---------------------------------
# The compare and contrast with the map.Reduce()function from the package.
# * Changed the data structure returned
# 	to a list (SIMPLIFY=FALSE)
#	This produces a list with length equal to the number
# 	of cells in the partition, a convenient
# 	option when we will need to pass each cell's 
# 	predictive model object to the reduce function.   
# * like mapReduce() a data object, D, is split into 
# 	list corresponding to parition cells and passed
# 	to the reduce function. 
# * Added a mechanism for passing parameters 
# 	to the reduce function, e.g. a list of model objects.  
# * The reduce parameters are specified by the user 
#	who can control which reduce parameters are 
# 	local (that is they may vary by partition cell) and which 
# 	apply globally across the partition.
#
# TODO
# ------
# * The mapReduce() function has options to parallelize 
# 	the jobs via papply or multicore. 
# 	Can this be done here? 
# * Its important to parallelize at this core level 
# 	so that we break up the data into the smallest possible 
# 	pieces (smaller memory footprint)  
#-------------------------------------------------------------------------------
vMR <- function(
	map.function,
	map.data, 
	map.parameters, 
	reduce.function, 
	reduce.data, 
	reduce.parameter.vectorize = NULL,
	...){ 
#-------------------------------------------------------------------------------
map <- map.function(map.data, map.parameters)
split.data <- split(reduce.data, map)
vectorize.args <- c("D")
if (!is.null(reduce.parameter.vectorize))
	vectorize.args <- c("D",reduce.parameter.vectorize)
vReduceFunction <- Vectorize(reduce.function, 
    						vectorize.args = vectorize.args,
    						SIMPLIFY=FALSE)
vMR_result <- vReduceFunction(D=split.data, ...)				
return(vMR_result)
#-------------------------------------------------------------------------------
}# END FUNCTION
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Exsemble of Partitions - Model Diagnostics - VI
#-------------------------------------------------------------------------------	
# For some prediction functions, like GAM, many 
# model diagnostics and statistical inferences are 
# already computed. So, we also create a 
# Ensemble Diagnostics step, as a user interface to 
# show how other model summaries and objects 
# can be computed and associated with the geo-indices. 
# e.g.  
#returned may include more than just the prediction, 
#eg. SE's and terms, etc. Model summaries and data objects
#like these can be collected by writing an 
#appropriate Ensemble-level diagnostic function 
#(very much like these) tweaked to handle the 
#specific data structures returned by vMR. 


