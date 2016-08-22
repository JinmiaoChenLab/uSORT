itn_getMatchingOriginalName<-function (new, original, ignore_case = FALSE) 
{
    if (ignore_case == FALSE) 
        return(new)
    org_sub <- original[itn_isElementIgnoreCase(original, new, 
                                                ignore_case)]
    if (length(org_sub) != length(new)) 
        stop("Error in itn_getOriginalName.")
    updated_new <- original[toupper(original) %in% toupper(new)]
    updated_new_mtx <- as.matrix(updated_new)
    rownames(updated_new_mtx) <- toupper(updated_new)
    updated_new <- as.character(updated_new_mtx[toupper(new), 
                                                ])
    return(updated_new)
}
