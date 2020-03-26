#' Check first column of matrix for names
#'
#' If values in first column are not integers or numerics, the first column will be changed to a row name.
#' @param df A data frame
#' @return A data frame with row names from the first column if the first column was not numeric or integer

check.name <- function(df){
	if(is.null(df)){
		return(df)
	}else{
		if (class(df[,1]) != "integer" & class(df[,1]) != "numeric"){
			df <- df[,-1]
		}
	}
	return(df)
	}
