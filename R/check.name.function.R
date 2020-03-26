#' Check first column of data frame for names
#'
#' If values in first column of data frame are not integers or numerics, the first column will be erased.
#' @param df A data frame
#' @return A data frame (with first column removed if values were not integers or numerics)

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
