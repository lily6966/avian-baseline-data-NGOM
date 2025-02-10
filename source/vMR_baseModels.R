



#-------------------------------------------------------------------------------
# baseModel = GAM
# Model Traing, Prediction & Diagnostic Functions
#------------------------------------------------------------------------------
# 
# The baseModel tasks (training, prediction, diagnostics) 
# takes place on a cell by cell basis. 
# For model training, this means that the data input 
# are a subset of the data for a single cell. 
# Similarly, the results from each baseModel task are passed 
# back to the MR-Ensemble wrapper on a cell-by-cell basis
# where they are assembled into Ensemble-level results. 
# For this reason output from model training and prediction must 
# adhere to the ensemble definition. 
#
# Input requirements for model training, prediction, and 
# diagnostics are spelled out below. The input requirements
# have been kept to a minium, and users can pass additional 
# input parameters. This means that users have a high 
# degree of flexibility writing new training, prediction, 
# functions. 
#
# 		For now parameters passed will be treated globally, i.e. 
# 		the same parameter value used for all cells in the partition. 
# 		vMR() may allow for "local" or cell-level or vectorized parameter values.
# 		The vectorized parameters need further development. Note,
# 		I do know that to use vectorized parameters, they MUST be 
# 		defined formally in the baseModel function defintiion 
# 		(not just ellipses) in order for the parameter vectorization to work. 
#
# From the coding perspective the difference between Training, 
# Prediction, and Diagnostics tasks is essentailly a difference
# between API's, i.e. differences in data structure for IO. 
# 
# TRAINING
#	I: D
#	O: list(model object)
# 	Ensemble Product: list of model object lists
# PREDICTION
# 	I: D, model.obj
# 	O: numeric (n x 1) prediction vector
# 	Ensemble Product: mean & sd prediction vectors
# DIAGNOSTICS
# 	I: D, model.obj
# 	O: list(user defined diagnostic object)
# 	Ensemble Product: list of user defined diagnostic objects
#
# TODO: 
# ---------
# I have not written the diagnostics module yet, but I 
# want to create VI & PD modules for GAM and GBM that 
# use baseModel routines to compute these local stats for me.
# The idea is that these local results will simply 
# be passes through an EnsebleDiagnostics function, like 
# model training. Then I can write a summary wrapper
# that will extract local VI's and another that will extract 
# local partial effects. 
# E.g. 
#
# I need to think about how to pass in weights for model training. 
# There are a few different ways to achieve this -
# simplest global parameter, then vectorized or local parameters. 
#
#-------------------------------------------------------------------------------
# GAM_train Documentation
#
# REQUIRED INPUT:
# ----------------
# 	D - a data frame, including both response and predictors. 
#	Remember that each baseModel will only see that data, i.e. 
# 	the subset of rows of D for a single cell at a time. 
#
# OPTIONAL INPUT:
# ----------------
# 	baseModel parameters can be passed in from user call to 
# 		ENSEMBLE level call. 
# 
# OUTPUT:
# -----------
# 	Essentially thre are two values passed back. 
# 	
#	1) In the case that the baseModel operation/task can not be completed
# 	the value list(NA) must be passed back to the ensemble wrapper.  
#	
#	2) When conditions are satistisfied to carry out baseModel task
# 	(training or prediction) a list pobject must be returned.
#	In practice most of the R model objects are package specific lists 
# 	or are easy to wrap up as lists.
# 	This is a very flexibile way to pass back data & results. 
#		??? have I tried otheres? or 
# 		must the returned object be a list? 
#	
# Discussion
# --------------  
# baseModelTraining is a good place to establish minimum data levels.
# added min.data.size parameter with a default value of 40. 
# In addition I have added a parameter to specify the formula. 
# But the rest of the gam parameters (family, etc.) are 
# specified by the user when they call the EnsembleModelTrain. 
#	
#-------------------------------------------------------------------------------
	GAM_train <- function(D,
						bm.formula,
						min.data.size = 40, 
						...){ 
	require(mgcv)
	d.gam <- list(NA)
	if (nrow(D) > min.data.size){
		d.gam <- gam(formula=as.formula(bm.formula),
						data=D, 
						...)			
		}
	return(d.gam)
	} 
		
#-------------------------------------------------------------------------------
# GAM_Predict 
#
# REQUIRED INPUT:
# ----------------
# 	D - a data frame, including both response and predictors. 
#		Remember that each baseModel will only see that data, i.e. 
# 		the subset of rows of D for a single cell at a time.
# 	model.obj -  
#
# OPTIONAL INPUT:
# --------------- 
# Added two parameters to insure that predict.gam 
# returns the vector of response scale predictions. 
# 	se.fit = F by default becuase data structure returned to 
# 				pred may not be a vector otherwise!!!
#	type="response". Note, by adding this to default 
# 		value it is necessary to also include the type=type \
#		on the call to pred. There are some subtlties about
# 		how parameter values are assigned and passed. 
# 			pred <- predict(model.obj, newdata=D, 
#				se.fit=se.fit, type=type, ... )	
#
# 	baseModel parameters can be passed in from user call to 
#	 ENSEMBLE level call. 
#
# OUTPUT:
# -----------
# 	numeric vector - will be reassembeled by ensemblepredict function.
#-------------------------------------------------------------------------------	
	GAM_predict <- function(D,model.obj,
			se.fit=F,
			type="response", ...){ 	
	require(mgcv)
	pred <- NULL
	# This condition generates a warning when there IS 
	# a gam object because its a list with length > 1
	# 4.21.10 - Now I check the first element of the list
	model.na <- length(model.obj)==1 & is.na(model.obj[[1]][1])		
	if ( !model.na ) 
		pred <- predict(model.obj, newdata=D, 
				se.fit=se.fit, type=type, ... )	
	return(pred)
	} 

# NOT TESTED
	GAM_VI <- function(D,model.obj,...){
	require(mgcv)
	model.na <- length(model.obj)==1 & is.na(model.obj[[1]][1])		
	if ( !model.na ){ 	
		ttt <- anova(model.obj)$chi.sq
		vi <- ttt/sum(ttt)*100		
		}
	return(list(vi=vi))
	} 


#-------------------------------------------------------------------------------
# GBM baseModel Train and Predict 
#
# NOTE: GBM will crash if ALL values for predictor are NA
#-------------------------------------------------------------------------------
	GBM_train <- function(D,
						bm.formula, 
						min.data.size = 40, 
						...){ 
	require(gbm)
	d.gbm <- list(NA)
	if (nrow(D) > min.data.size){
		d.gbm <- gbm(formula=as.formula(bm.formula),
						data=D, 
						...)			
		}
	return(d.gbm)
	} 
	# -----------------------------------------------
	# Its important for the user to pass n.trees!!!!
	GBM_predict <- function(D,model.obj,...){ 	
	require(gbm)
	pred <- NULL
	# This condition generates a warning when there IS 
	# a gam object because its a list with length > 1
	# 4.21.10 - Now I check the first element of the list
	model.na <- length(model.obj)==1 & is.na(model.obj[[1]][1])		
	if ( !model.na ) 
		pred <- predict(model.obj, newdata=D, ...)	
	return(pred)
	} 

#-------------------------------------------------------------------------------
# GBM baseModel Training
# Modified for Predictors that may be uniquely NA 
# 7.11.11
# 
# This version checks all predictors for NA's 
# and removes those predictors that are uniquely NA
# to avoid the error - GBM will crash if ALL values for 
# predictor are NA.
#
# Note that the previous basemodel prediction function 
# works fine with these model objects. predict.gbm 
# will ignore extra variables in the prediction design. 
#
# Input 
# --------
# * The response variable MUST be named "y"
# * There is NO argument bm.formula
#-------------------------------------------------------------------------------
	GBM_NAremove_train <- function(D,
						min.data.size = 40, 
						...){ 
	require(gbm)
	d.gbm <- list(NA)
	if (nrow(D) > min.data.size){
		D <- D[ , apply(is.na(D), 2, sum) != nrow(D)]  
		non.na.pred.names <- names(D)[ names(D) != "y"]
		model.formula <- paste(non.na.pred.names, collapse="+")
		model.formula <- paste("y ~ ", model.formula)
		d.gbm <- gbm(formula=as.formula(model.formula),
						data=D, 
						...)			
		}
	return(d.gbm)
	} 


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# stem.gbm.predictor.importance 
# 7.14.11
#
# For each partition we create a data frame containing 
# the stixel location and the predictor importance (PI)s
# from the base modle in that location.  
#
# It will be important to assume that each PI vector may have 
# missing predictors. We add in a step that checks for 
# missing predictors and packs the results matrix against 
# as master set of predictors.
#
# INPUT:
# -------------
# 	stem.object - the model object
# 	partition.number = 1, currently this function computes
# 		PI's one partition at a time
#	n.trees - the gbm parameter
#	predictor.names - vector of character strings that is 
# 		the master predictor set. 
#
# OUTPUT:
# ------------------ 
# The columns of the data frame are:
#		x.centroid
#       y.centroid
#		t.centroid
# Followed by one column for each of the predictors used in the 
# design matrix. Each row contains the PI for a single base model 
# supported by a single stixel. The row order should be 
# the same as the parition order. 
#
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
stem.gbm.predictor.importance <- function(
	stem.object,
 	partition.number = 1,
	n.trees, 
	predictor.names) {
# -------------------------------
# ------------------------------- 
pn <- partition.number 
pi.data <- matrix(NA, 
	length(stem.object$model[[pn]]), 
	length(predictor.names) )
pi.locs <- matrix(NA,
	length(stem.object$model[[pn]]), 
	ncol(stem.object$partition[[pn]]$min))
pi.data <- as.data.frame(pi.data)
names(pi.data) <- c(predictor.names )
pi.locs <- as.data.frame(pi.locs)
names(pi.locs) <- paste(names(stem.object$partition[[pn]]$min),
			".centroid",sep="")	 
# Loop over all base-models
for (iii in 1:length(stem.object$model[[1]]) ){
	#iii <- 102
	mn <- iii 	# model object number			
	partition.cell.number <- as.numeric(names(stem.object$model[[pn]])[mn]) 
	base.model.obj <- stem.object$model[[pn]][mn][[1]]
	model.na <- length(base.model.obj)==1 & 
					is.na(base.model.obj[[1]][1])		
	#if (model.na) 
	#	cat(" base.model from cell ", partition.cell.number, " is NA \n") 
	if (!model.na){
		# Compute Partition Cell Centroid	
			pi.locs[iii,] <- 
				(stem.object$partition[[pn]]$min[partition.cell.number,] +
				stem.object$partition[[pn]]$max[partition.cell.number,] )/2
		# Compute Predictor Importance 
			ttt.df <- summary.gbm(base.model.obj,
			        n.trees = n.trees,
			        plotit = FALSE,
			        order = FALSE,
			        method = relative.influence,
			        normalize = FALSE)       
		# Repack predictor importance
		ttt.ind <- match(predictor.names, as.character(ttt.df$var)) 
		pi.data[iii,] <-  ttt.df$rel.inf[ttt.ind]         
		# -------------------------------
	} # end if (!model.na){		
} # end iii	
# ---------------------------------------------	
return(data.frame(pi.locs,pi.data))
} #end function
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# RPART baseModel Train and Predict 
# 
# parameter "method" is used for response family
# rpart.control(cp=0.01,xval=0,minbucket=10))
#-------------------------------------------------------------------------------
	RPART_train <- function(D,
						bm.formula, 
						min.data.size = 20, 
						...){ 
	require(rpart)
	d.rpart <- list(NA)
	if (nrow(D) > min.data.size){
		d.rpart <- rpart(formula=as.formula(bm.formula),
						data=D, 
						...)			
		}
	return(d.rpart)
	} 

	RPART_predict <- function(D,model.obj,...){ 	
	require(rpart)
	pred <- NULL
	# This condition generates a warning when there IS 
	# a gam object because its a list with length > 1
	# 4.21.10 - Now I check the first element of the list
	model.na <- length(model.obj)==1 & is.na(model.obj[[1]][1])		
	if ( !model.na ){ 
		# -------------------------------------------------
		# rpart Predict
		# ------------------
		# Note: the predict method for rpart will
		# not make a prediction on a root node only.
		# That is, if the rpart model is only a root node
		# then I need to skip the call to predict - will generate an error
		# and manually assign all pts the same prediction.
		# --------------------
		# 4.9.09
		# Separate out case when bernoulli because
		# it requires slightly different call to deal with two classes.
		# ---------------------
		# 1.31.11
		# I am using the same "trick" as the stem.library.b.R
		# subtracting 1- model.obje$frame$yval
		# I do not think that this will work with other 
		# codings of the response variable levels, 
		# or other types of response variable. 
		# -------------------------------------------------
		#NROW(f.rpart$frame)
		if (NROW(model.obj$frame) == 1 & 
			model.obj$method == "class") {
				# This is the function that may rely on
				# the specific coding used for binary class
				# classification. We are using FACTORS
				# "TRUE" and "FALSE". The alphanumeric values 
				# are FALSE = 1 and TRUE = 2. 
			    pred <- rep((model.obj$frame$yval-1), NROW(D))
		} #end class& stump predictions section
		if (NROW(model.obj$frame) > 1 | 
			model.obj$method != "class") {	
			# just use the regular prediction 
			pred <- predict(model.obj, newdata=D, ...)	
			# Assume only two classes and send back second.
			# This means that "success" must be coded with 
			# "largest" class. Largest probably means 
			# the largest in terms of as.numeric() etc.  
			if (model.obj$method == "class") pred <- pred[,2]
		}#end class& stump predictions section
	}	
	return(pred)
	} 


#--------------------------------------------------------------------------
# Support Function
#--------------------------------------------------------------------------
	support_train <- function(D, 
						min.data.size = 0){ 
	d.mod <- list(NA)
	if (nrow(D) > min.data.size){
		d.mod <- TRUE			
		}
	return(d.mod)
	} 
	# --------------------------------------------
	support_predict <- function(D,model.obj){ 	
	pred <- NULL
	model.na <- length(model.obj)==1 & is.na(model.obj[[1]][1])		
	if ( !model.na ) 
		pred <- 1
	return(pred)
	} 



