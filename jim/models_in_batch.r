# An alternative way to build regression models in a batch

# ys is the matrix of multiple dependant variables
# xs is the matrix of multiple independant variables

ys.list <- lapply(1:ncol(ys), function(i) ys[, i]); # turn ys into a list;

d <- data.frame(xs, one_y = ys.list[[1]]); # create the data.frame for modeling

# Create the formula
e <- paste(colnames(xs), collapse=' + '); 
f <- formula(paste('one_y', e, sep=' ~ '));

result <- parallel::mclapply(ys.list, function(y) { # Run parallel computing
  
  d$one_y <- y; # replace the dependent variable with a new one
  
  # create model using glm or speedglm
  model <- glm(f, data=d, family=binomial);
  
  # retrieve statistics from the model, such as p value of Wald test
  c(); # return the statistics
  
}, mc.core=4); 

