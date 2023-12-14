matrix.split.train.test <- function(mask, testp=0.2, by.row=FALSE, seed=2023 ) {
   #' returns a new mask similar to mask with a new train-test sets.
   #' testp: proportion of test out of the nonzero cells in mask
   #' new_mask equals 1 for the training and for when mask is set to 0. 
   #' It equals 0 for the new test split.
   
   set.seed(seed)
   n_rows <- dim(mask)[1]
   n_cols <- dim(mask)[2]
   # Create a data frame of all matrix indices
   # we only consider non-missing data (ie, with mask_ij=1)
   indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)[mask==1,]
   # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
   indices <-  indices[sample(1:nrow(indices)), ]
   row.names(indices) <- NULL
   
   test.idx = sample(1:nrow(indices),size = nrow(indices)*testp, replace = FALSE)
   test.indices = indices[test.idx,]
   
   new_mask = matrix(1, nrow=n_rows, ncol = n_cols)
   new_mask[as.matrix(test.indices[, c("row", "col")])] <- 0
   
   return(new_mask)
}

