############################################################################################
##                              Library Installation and Loaded                           ##
############################################################################################

packagesLoad <- function(pack){
  packages<-foreach(i=1:(dim(pack)[1]), .combine='c') %dopar% paste(pack[i,(dim(pack)[2])])
  return(packages)
}

packagesInstall <- function(packages){
    new.pkg <- packages[!(packages %in% installed.packages()[, "Package"])]
    if(length(new.pkg)){ 
        install.packages(new.pkg, dependencies = TRUE)
    }
    sapply(packages, require, character.only = TRUE)
}

############################################################################################