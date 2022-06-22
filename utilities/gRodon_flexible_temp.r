#### Modified from code contributed by Jake Weissman ####

## Don't forget to set your working directory!

## User needs to set these variables

prefix <- 'test.archaea'      # prefix for your analysis files
paprica.path <- '~/paprica'  # assumes paprica is in your home directory, change if necessary

## Functions

boxcoxTransform <- function(x, lambda, back_transform = F) {
  #This function needed to back-transform gRodon's raw model estimates out of
  #box-cox tranformed space
  if (back_transform == TRUE) {
    return((x*lambda +1)^(1/lambda))
  } else {
    return((((x^lambda) - 1) / lambda))
  }
}

predictFromModel <- function(x,model,lambda_model){
  #Predicts from model and reverses box-cox transform
  return(boxcoxTransform(predict.lm(model,x),
                         lambda_model,
                         back_transform = T))
}

## Load necessary gRodon models

load(paste0(paprica.path, '/models/', 'gRodon_models.rda'))

## Load edge_data.csv file and modify as necessary

edge.data <- read.csv(paste(prefix, 'edge_data.csv', sep = '.'))
edge.data.grodon <- edge.data[,c('gRodon.CUBHE','gRodon.ConsistencyHE','gRodon.CPB')]
colnames(edge.data.grodon) <- c('CUBHE', 'ConsistencyHE', 'CPB')

## Add a column with optimal growth temps.  This example code assumes a common growth optimum for the sample,
## but this could vary by edge.

edge.data.grodon$OGT <- rep(12, dim(edge.data.grodon)[1])

predictFromModel(edge.data.grodon[,c('CUBHE','ConsistencyHE','CPB', 'OGT')],
                 gRodon_model_temp_madin,
                 lambda_milc_madin)
