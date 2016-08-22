itn_isSubSetIgnoreCase<- function (subset, set, ignore_case = FALSE) 
{
    subset <- as.vector(subset)
    set <- as.vector(set)
    if (length(set) == 0 || length(subset) == 0) {
        return(FALSE)
    }
    if (ignore_case == TRUE) {
        subset <- tolower(subset)
        set <- tolower(set)
    }
    if (length((intersect(subset, set))) == length(subset)) {
        return(TRUE)
    }
    else {
        return(FALSE)
    }
}