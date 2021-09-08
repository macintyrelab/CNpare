######################################
## FUNCTIONS FOR COMPARING PROFILES ##
######################################


#' @title Perform a pearson correlation test
#' @description This function calculates the Pearson's r
#' @name pair_pearson
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of model
#' @param ccle.name name of model
#'
#' @return Pearson's r
#' @export
#' @examples

pair_pearson<-function(cell, ccle, ccle.name){
    cor <- cor(cell, ccle, use = "na.or.complete", method = "pearson")
    cor <- as.data.frame(cbind(id=ccle.name,
                              cor.coef=cor))
    rownames(cor) <- NULL
    return(cor)
}


#' @title Calculation of manhattan distance
#' @description This function calculates the Manhattan distance
#' @name pair_manhattan
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of model
#' @param ccle.name name of model
#'
#' @return Manhattan distance
#' @export
#' @examples

pair_manhattan<-function(cell, ccle, ccle.name){
    dist <- cbind(id=ccle.name,
                 distance=mean(abs(cell - ccle),na.rm=TRUE))
    rownames(dist) <- NULL
    return(dist)
}


#' @title Calculation of euclidean distance
#' @description This function calculates the Euclidean distance
#' @name pair_euclidean
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of model
#' @param ccle.name name of model
#'
#' @return Euclidean distance
#' @export
#' @examples

pair_euclidean<-function(cell, ccle, ccle.name){
    dist <- cbind(id=ccle.name,
                  distance=sqrt(sum((cell - ccle)^2,na.rm=TRUE)))

    rownames(dist) <- NULL
    return(dist)
}


#' @title Calculation of cosine similarity
#' @description This function calculates the Cosine similarity
#' @name pair_cosine
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of model
#' @param ccle.name name of model
#'
#' @return Cosine similarity
#' @export
#' @examples

pair_cosine<-function(cell, ccle, ccle.name){
    cos_sim <- (sum(cell*ccle,na.rm=TRUE))/
        (sqrt(sum(cell^2,na.rm = TRUE))*sqrt(sum(ccle^2,na.rm = TRUE)))
    dist <- cbind(id=ccle.name,
                  cos_sim=cos_sim,
                  distance=1-cos_sim)

    rownames(dist) <- NULL
    return(dist)
}


#' @title Get similarity metrics
#' @description This function get similarity metrics for comparisons
#' @name getSimilarities
#'
#' @param dat1 matrix with bin-level copy numbers per sample to test
#' @param dat2 matrix with bin-level copy numbers per model to compare
#' @param method method to use for testing similarity. Default is all
#' @return dataframe with similarity metrics for all comparisons
#' @export
#' @examples

getSimilarities<-function(dat1, dat2, method="all"){
    #get matrix with cell lines copy-number
    cn_matrix<-getInputmatrix(dat1)

    #get list with correlation matrix
    out   <- list()
    samps <- colnames(cn_matrix)
    for (i in samps)
    {
        cell_cn  <- cn_matrix[,which(colnames(cn_matrix)==i)]

        if (method == "all" | method == "pearson"){
            corr <- lapply(1:ncol(dat2),
                           function(x) pair_pearson(cell=cell_cn,
                                                         ccle=dat2[,x],
                                                         ccle.name=colnames(dat2)[x]))
            corr <- do.call(rbind,corr)
            out[['pearson']][[i]] <- corr

        }
        if (method == "all" | method == "manhattan"){
            man <- lapply(1:ncol(dat2),
                           function(x) pair_manhattan(cell=cell_cn,
                                                      ccle=dat2[,x],
                                                      ccle.name=colnames(dat2)[x]))
            man <- do.call(rbind,man)
            out[['manhattan']][[i]] <- as.data.frame(man)
        }
        if (method == "all" | method == "euclidean"){
            eu <- lapply(1:ncol(dat2),
                          function(x) pair_euclidean(cell=cell_cn,
                                                     ccle=dat2[,x],
                                                     ccle.name=colnames(dat2)[x]))
            eu <- do.call(rbind,eu)
            out[['euclidean']][[i]] <- as.data.frame(eu)
        }
        if (method == "all" | method == "cosine"){
            cos <- lapply(1:ncol(dat2),
                         function(x) pair_cosine(cell=cell_cn,
                                                    ccle=dat2[,x],
                                                    ccle.name=colnames(dat2)[x]))
            cos <- do.call(rbind,cos)
            out[['cosine']][[i]] <- as.data.frame(cos)
        }
    }

    pearson<-data.table::rbindlist(out[['pearson']], idcol=TRUE)
    pearson$cor.coef<-as.numeric(pearson$cor.coef)
    manhattan<-data.table::rbindlist(out[['manhattan']], idcol=TRUE)
    manhattan$distance<-as.numeric(manhattan$distance)
    euclidean<-data.table::rbindlist(out[['euclidean']], idcol=TRUE)
    euclidean$distance<-as.numeric(euclidean$distance)
    cosine<-data.table::rbindlist(out[['cosine']], idcol=TRUE)
    cosine$cos_sim<-as.numeric(cosine$cos_sim)
    cosine$distance<-as.numeric(cosine$distance)

    #join all metrics and return a data.frame
    if (method=="all"){
        measures<-cbind(pearson,manhattan[,3],euclidean[,3],cosine[,3])
        colnames(measures)<-c("fileid", "id","r","manhattan","euclidean","cos_sim")
        return(measures)
    }

    #or return only the metric desired
    else if (method=="pearson"){
        return(pearson)
    }
    else if (method=="manhattan"){
        return(manhattan)
    }
    else if (method=="cosine"){
        return(euclidean)
    }
    else if (method=="cosine"){
        return(cosine)
    }
}


#' @title Get top hits
#' @description This function get the top hit for each comparison per sample
#' @name getTopHit
#'
#' @param samples vector with name of samples
#' @param measure dataframe with similarity metrics
#' @param method similarity metric to use for ordering . Default is all
#' @return dataframe with similarity metrics for all comparisons
#' @export
#' @examples

getTopHit<-function(samples, measure, method="all"){
    colnames(measure)[1]<-"cellid"
    tops<-list()
    if (method=="all" | method=="pearson"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(-m$r),]
            top<-m[1,c(1:2)]
            tops[['pearson']][[i]]<-top
        }
    }
    if (method=="all" | method=="cosine"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(-m$cos_sim),]
            top<-m[1,c(1:2)]
            tops[['cosine']][[i]]<-top
        }
    }
    if (method=="all" | method=="manhattan"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(m$manhattan),]
            top<-m[1,c(1:2)]
            tops[['manhattan']][[i]]<-top
        }
    }
    if (method=="all" | method=="euclidean"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(m$euclidean),]
            top<-m[1,c(1:2)]
            tops[['euclidean']][[i]]<-top
        }
    }

    pearson<-data.table::rbindlist(tops[['pearson']])
    manhattan<-data.table::rbindlist(tops[['manhattan']])
    euclidean<-data.table::rbindlist(tops[['euclidean']])
    cosine<-data.table::rbindlist(tops[['cosine']])

    #join all top hits and return a data.frame
    if (method=="all"){
        tophits<-cbind(pearson[,c(1,2)],cosine[,2],manhattan[,2],euclidean[,2])
        colnames(tophits)<-c("id", "topPearson","topCosine", "topManhattan","topEuclidean")
        return(tophits)
    }

    #or return only the metric desired
    else if (method=="pearson"){
        return(pearson)
    }
    else if (method=="manhattan"){
        return(manhattan)
    }
    else if (method=="cosine"){
        return(euclidean)
    }
    else if (method=="cosine"){
        return(cosine)
    }
}



#' @title Get input data as matrix
#' @description This function prepare data for comparisons
#' @name getInputmatrix
#'
#' @param dat bin-level copy numbers
#' @return matrix with copy numbers of all samples
#' @export
#' @examples

getInputmatrix<-function(dat){
    #prepare matrix with cell lines copy-number
    if(class(dat) == "QDNAseqCopyNumbers") {
        cn_matrix <- dat@assayData$copynumber
    } else {
        cn_matrix <- dat
    }
    return(cn_matrix)
}
