#included are the checkpoints, warnings and stop functions
#performs the mvGSEA function in parallel

edited_wrapper_mvGSEA <- function(ys,X, v, conf=NULL, mc=4, prediction_method=c("sensitivity", "custom"),prediction_cutoff=0.005)
{
  
  #check that v and conf don't overlap
  if (!is.null(conf)){
  if (length(Reduce(intersect, list(v, conf))) != 0){
    stop("v and conf inputs must be distinct")}
  }

  
  
#------------------------------------------------------------
#if prediction_cutoff is not a decimal, then set it to the default
  if(prediction_cutoff <= 0 | prediction_cutoff >= 1 ){
    prediction_cutoff=0.5
  }
  
#------------------------------------------------------------
  #if mc is some strange number, set it to the default number of cores
 if (!is.element(mc,seq(from = 2, to = 64, by = 2))){
    mc=4
  }

  #------------------------------------------------
  

  if (length(rownames(X))==0 | length(rownames(ys))==0) {
    #if there are no row names, AND if the number of rows of ys and X aren't equal, put in a warning to specify row names
    #if there aren't row names, but if X and ys contain equal number of rows, then otherwise assume the rows match up
    if(length(row(ys)) != length(row(X))){
    stop("Must specify row names")
    } else{ X_new_ <- X
    ys_new_ <- ys
    }
  } else{ #if there ARE row names
  #if row names of y and x don't match, take out the ones that don't overlap, and then order the ones that match
  common <- Reduce(intersect, list(rownames(X), rownames(ys))) #the common ones
  
 #--------------------------------------
   #if there aren't any COMMON row names
  if (length(common)<3)
    stop("There are less than 3 row names. At least 3 common rows are required.")
  if (length(common) < 20){
  warning("There less than 20 common genes. Need more for meaningful analysis.")}
  
  #-------------------------------------------------
   #take out the ones that AREN'T the common ones
   #for X, taking out the non-overlapping row names-----
  X_indices <- pmatch(common, rownames(X)) #the indices that are common in the X matrix
  X_new <- X[X_indices,,drop=FALSE] #using the indices to retrieve the common index VALUES(row names)
  X_new_ <- X_new[gtools::mixedsort(rownames(X_new), decreasing = FALSE), , drop=FALSE] #take these VALUES(row names) and order them
  #for ys, taking out the non-overlapping row names-----
  ys_indices <- pmatch(common,rownames(ys)) #the indices that are common in the X matrix
  ys_new <- ys[ys_indices,,drop=FALSE] #using the indices to retrieve the common index VALUES(row names)
  ys_new_ <- ys_new[gtools::mixedsort(rownames(ys_new), decreasing = FALSE), , drop=FALSE] #take these VALUES(row names) and order them
  }
  
  #------------------------------------------------------
  
  
  #if there's NA, set it to 0    
  ys_new_[is.na(ys_new_)]<-0 #assume missing value is 0
  
  #if ys values aren't binary
  ys_new_bin <- ys_new_[ys_new_ != 0 & ys_new_ != 1]
  ys_new_bin
  if (length(ys_new_bin) != 0) {
    stop("ys must contain only binary values")
  }
  
  #-----------------------------------------
  
  #check that v and conf names are in X (columns)
  check_v_in_X <- v %in% colnames(X_new_)
  

  
  check_conf_in_X <- c()
  if (!is.null(conf)){
    check_conf_in_X <- conf %in% colnames(X_new_)
  }else check_conf_in_X <- TRUE
  
  if(all(check_v_in_X==TRUE)==FALSE | all(check_conf_in_X==TRUE)==FALSE)
    stop("v or conf columns aren't found in the X matrix")
  
  #extracting the v columns
  selected_v_indices <- pmatch(v, dimnames(X_new_)[[2]])
 selected_v_columns <- X_new_[,selected_v_indices, drop=FALSE]
  
  #extracing the confounder columns
  confounders_indices <- c()
  confounders_columns <- c()
  if( !is.null(conf)){ #if there ARE confounders
  confounders_indices <- pmatch(conf, dimnames(X_new_)[[2]]) #find which columns are confounders
    confounders_columns <- X_new_[,confounders_indices, drop=FALSE]}#confounder column(s)
  
 #----------------------------------------
  
ys.list <- lapply(1:ncol(ys_new_), function(i) ys_new_[,i])

if (is.null(confounders_columns)){
d <- data.frame(selected_v_columns, one_y = ys.list[[1]])
if (length(colnames(selected_v_columns)) > 1){
v_prod <- paste(colnames(selected_v_columns), collapse=' * ')
} else{
  v_prod <- colnames(selected_v_columns)
}
f <- formula(paste('one_y', v_prod, sep=' ~ '))
partial <- NULL
} 
if (!is.null(confounders_columns)){ #if there ARE CONFOUNDERS
 d <- data.frame(selected_v_columns, confounders_columns, one_y = ys.list[[1]])
  v_prod <- paste(colnames(selected_v_columns), collapse=' * ')
  if (length(colnames(confounders_columns))>1){
    conf_prod <- paste(colnames(confounders_columns), collapse= '*') 
  } else{ #if there is only one confounder
    conf_prod <- colnames(confounders_columns)
  }
  f <- formula(paste('one_y', paste(v_prod,conf_prod, sep='+'), sep=' ~ '))
  partial <- formula(paste('one_y', conf_prod, sep='~'))
}



result <- parallel::mclapply(ys.list, function(y){
    
  d$one_y <- y
    mvGSEA(y,X=X_new_,v_col=selected_v_columns,conf_col=confounders_columns, prediction_method=prediction_method,prediction_cutoff=prediction_cutoff, f=f, d=d, partial=partial)
    
    
  }, mc.cores=mc);


#----------------------------------------- 
 output_table <- Reduce(rbind, result)
  row.names(output_table) <- colnames(ys)
  #return a table with each column in ys as a row
 output_table 


}
