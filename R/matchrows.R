#' @title matchrows
#' @description groups rows into a new column, a number of strings is used to find matching rows and marks them with a tag
#' #'
#' @param df dataframe
#' @param colname name of the column to be grouped, exe: colname="Species"
#' @param newcol name of the new column, should be a different name than the existing columns, default set to "groups"
#' @param string vector with the strings to match, matching is made by grepl
#' @param tag name of the matching rows that will appear in the newcol, default set to "grouped"
#' @examples
#' \dontrun{
#' quiet(cat("test"))
#' quiet(warning("test"))
#' quiet(warning("test"), all=T)
#' }
#' # This is a function that suppresses log info
#' @export
matchrows<-function(df, colname, newcol="groups", string=c(), tag="grouped" ){

  for (n in 1:length(string)){
    for (i in 1:nrow(df)){
      if(grepl(string[n], df[[colname]][i], fixed=TRUE)==TRUE){#Finds rows that match the substrings

        df[[newcol]][i]<-tag #string to substitute

      }else{}
    }#end of for loop i
  }#end of for loop n

  return(df)
}#end of function
