

#############################################################################################
#############################     train.test    #############################################
#############################################################################################

train.test <- function(data, perc=0.75, seed=100)
            {
            set.seed(seed)
            smp_size <- floor(perc * nrow(data))
            train_ind <- sample(seq_len(nrow(data)), size = smp_size, replace=FALSE)
            train <- data[train_ind, ]
            test  <- data[-train_ind, ]
            list  <- list()
            list  <- list(train=train,test=test)
            return(list)
 }













