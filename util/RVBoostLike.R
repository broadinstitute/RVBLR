#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(tidyverse))


parser = ArgumentParser()
parser$add_argument("--feature_matrix", required=TRUE, nargs=1, help="feature matrix")
parser$add_argument("--attributes", default="DJ,ReadPosRankSum,QD,FS,ED,PctExtPos", nargs=1, help="columns to boost on, comma-delimited")
parser$add_argument("--boosting_score_threshold", default=0.05, nargs=1, help="ecdf cutoff for defining true variants")
parser$add_argument("--output", required=TRUE, nargs=1, help="output matrix containing boosting results")

args = parser$parse_args()


###arguments
input_matrix = args$feature_matrix
output = args$output
attributes = strsplit(args$attributes, ",")[[1]]

if (! 'RS' %in% attributes) {
    attributes = c(attributes, 'RS')  # required
}

if (! 'IND' %in% attributes) {
    attributes = c(attributes, 'IND')  # required
}

boosting_score_threshold = args$boosting_score_threshold
result_table = args$output

message("Using attribute list: ", paste(attributes, collapse=","))

library(gbm)

data = read.table(input_matrix, header=T, stringsAsFactors = FALSE, sep="\t")


if (! all(attributes %in% colnames(data))) {
    missing_atts = attributes[ ! attributes %in% colnames(data) ]
    stop(paste("Error, missing attributes in data matrix:", missing_atts, sep=","))
}



# remove rna-editing sites
if ("RNAEDIT" %in% colnames(data)) {

    data = data %>% select(-RNAEDIT)
    message("-removing RNAEDIT sites")
}


# pull just attribute columns
data = data[, attributes, drop=F] # restrict to what we want to analyze here and reorder columns


# get RS and chrpos indicators

RS = data %>% pull(RS)

chrpos = data %>%  pull(IND)

data = data %>% select(-c(IND, RS))


## reset attributes sans RS
attributes = colnames(data)


# Set values with NA to the median value
for(j in 1:ncol(data)){
    is_na <- which(is.na(data[ ,j]))

    if(length(is_na) > 0){
        # get the median of the non NA values
        not_na <- which(!is.na(data[ ,j]))
        median_value <- median(data[not_na,j])
        # set the Na's to the median
        data[is_na,j] <- median_value
    }
}
############################
## Run adaboost

message("Running adaboost - rvboost-style")

NUMTREES = 2e4  #rvboost defaults

data = data.matrix(data)

res = gbm.fit(x=data,
              y=RS,
              n.trees=NUMTREES,
              interaction.depth=2,
              distribution="adaboost",
              verbose=FALSE)


## as per rvboost:
## convert it to the 0-1 scale since the adaboost method gives the predictions on logit scale.
## http://stats.stackexchange.com/questions/37497/how-to-use-r-gbm-with-distribution-adaboost

fitted_values <-  plogis(res$fit)

ecdf_func <- ecdf(fitted_values[which(RS==1)])

# apply all scores to ecdf_func
fitted_value_scores = ecdf_func(fitted_values)

result_table = data.frame('chr:pos'=chrpos, RS=RS, RVBfitval=fitted_values, RVBscore=fitted_value_scores, boosted = (fitted_value_scores >= boosting_score_threshold), check.names=FALSE )

write.table(result_table, file=output, quote=F, row.names=F, sep="\t")
