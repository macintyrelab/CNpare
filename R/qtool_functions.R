######################################
## FUNCTIONS FOR COMPARING PROFILES ##
######################################


#' @title Perform a pearson correlation test
#' @description This function calculates the Pearson's r
#' @name pair_pearson
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of the cell model
#' @param ccle.name name of the cell model
#' @param pvalue a logical value to indicate if the empirical pvalue is calculated or not
#'
#' @return Pearson's r
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' mod_cell<-as.matrix(ccle_cn[,2])
#' r<-pair_pearson(exp_cell,mod_cell,unique(cells_segcn$sample)[2],pvalue=FALSE)


pair_pearson<-function(cell, ccle, ccle.name, pvalue){
    r<-cor(cell, ccle, use = "na.or.complete", method = "pearson")
    cor<-as.data.frame(cbind(id=ccle.name, r=r))

    if(pvalue==TRUE){
        ind<-Position(function(x) x < r, all_pearson)-1
        cor<-as.data.frame(cbind(id=ccle.name, r=r,r.pval=ind/length(all_pearson)))
    }
    rownames(cor) <- NULL
    return(cor)
}


#' @title Calculation of manhattan distance
#' @description This function calculates the Manhattan distance
#' @name pair_manhattan
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of the cell model
#' @param ccle.name name of the cell model
#' @param pvalue a logical value to indicate if the empirical pvalue is calculated or not
#'
#' @return Manhattan distance
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' mod_cell<-as.matrix(ccle_cn[,2])
#' m<-pair_manhattan(exp_cell,mod_cell,unique(cells_segcn$sample)[2],pvalue=FALSE)

pair_manhattan<-function(cell, ccle, ccle.name, pvalue){
    manhattan=mean(abs(cell - ccle),na.rm=TRUE)
    dist<-cbind(id=ccle.name, manhattan=manhattan)

    if(pvalue==TRUE){
        ind<-Position(function(x) x > manhattan, all_manhattan)-1
        dist<-cbind(id=ccle.name, manhattan=manhattan,m.pval=ind/length(all_manhattan))
    }
    rownames(dist) <- NULL
    return(dist)
}


#' @title Calculation of euclidean distance
#' @description This function calculates the Euclidean distance
#' @name pair_euclidean
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of the cell model
#' @param ccle.name name of the cell model
#' @param pvalue a logical value to indicate if the empirical pvalue is calculated or not
#'
#' @return Euclidean distance
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' mod_cell<-as.matrix(ccle_cn[,2])
#' e<-pair_euclidean(exp_cell,mod_cell,unique(cells_segcn$sample)[2],pvalue=FALSE)

pair_euclidean<-function(cell, ccle, ccle.name, pvalue){
    euclidean=sqrt(sum((cell - ccle)^2,na.rm=TRUE))
    dist<-cbind(id=ccle.name, euclidean=euclidean)
    if(pvalue==TRUE){
        ind<-Position(function(x) x > euclidean, all_euclidean)-1
        dist<-cbind(id=ccle.name, euclidean=euclidean, e.pval=ind/length(all_euclidean))
    }
    rownames(dist) <- NULL
    return(dist)
}


#' @title Calculation of cosine similarity
#' @description This function calculates the Cosine similarity
#' @name pair_cosine
#'
#' @param cell vector with bin-level copy number values of one sample
#' @param ccle vector with bin-level copy number values of the cell model
#' @param ccle.name name of the cell model
#' @param pvalue a logical value to indicate if the empirical pvalue is calculated or not
#'
#' @return Cosine similarity
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' mod_cell<-as.matrix(ccle_cn[,2])
#' c<-pair_cosine(exp_cell,mod_cell,unique(cells_segcn$sample)[2],pvalue=FALSE)

pair_cosine<-function(cell, ccle, ccle.name, pvalue){
    cos_sim<-(sum(cell*ccle,na.rm=TRUE))/
        (sqrt(sum(cell^2,na.rm = TRUE))*sqrt(sum(ccle^2,na.rm = TRUE)))
    dist<-cbind(id=ccle.name, cos_sim=cos_sim)

    if(pvalue==TRUE){
        ind<-Position(function(x) x < cos_sim, all_cosine)-1
        dist<-cbind(id=ccle.name, cos_sim=cos_sim, c.pval=ind/length(all_cosine))
    }
    rownames(dist) <- NULL
    return(dist)
}


#' @title Get similarity metrics
#' @description This function get similarity metrics for comparisons
#' @name getSimilarities
#'
#' @param dat1 matrix with bin-level copy numbers per sample to test
#' @param dat2 matrix with bin-level copy numbers per model to compare
#' @param method method to use for testing similarity. Default is *all* methods.
#' @param pvalue a logical value to indicate if the empirical pvalue is calculated or not
#'
#' @return dataframe with similarity metrics for all comparisons
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr) getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' colnames(exp_cell)<-unique(cells_segcn$sample)[1]
#' measures<-getSimilarities(dat1=exp_cell,dat2=ccle_cn,method="pearson",pvalue=FALSE)

getSimilarities<-function(dat1, dat2, method="all", pvalue=FALSE){
    samples <- colnames(dat1)
    out <- lapply(samples, function(s) getSimilarities.sample(s, dat1, dat2, method, pvalue))
    out <- as.data.frame(do.call(rbind,out))
    ncol <- seq(ncol(out),from=3,by=1)
    cols <- colnames(out)[ncol]
    out[cols] <- as.data.frame(sapply(out[cols], as.numeric))
    return(out)
}

#' @title Get similarity metrics per sample
#' @description This function get similarity metrics for comparisons per sample
#' @name getSimilarities.sample
#'
#' @param s sample name
#' @param dat1 matrix with bin-level copy numbers per sample to test
#' @param dat2 matrix with bin-level copy numbers per model to compare
#' @param method method to use for testing similarity. Default is *all* methods.
#' @param pvalue a logical value to indicate if the empirical pvalue is calculated or not
#'
#' @return dataframe with similarity metrics for all comparisons of the sample
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins,data=cells_segcn,samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' colnames(exp_cell)<-unique(cells_segcn$sample)[1]
#' measures<-getSimilarities.sample(s=colnames(exp_cell),dat1=exp_cell,
#'     dat2=ccle_cn,method="pearson",pvalue=FALSE)

getSimilarities.sample<-function(s, dat1, dat2, method, pvalue){
    #get matrix with cell lines copy-number
    cn_matrix <- getInputmatrix(dat1)

    #get list with correlation matrix
    out <- list()
    cell_cn <- cn_matrix[,which(colnames(cn_matrix)==s)]
    if (method == "all" | method == "pearson"){
        corr <- lapply(seq_len(ncol(dat2)),
                       function(x) pair_pearson(cell=cell_cn,
                                                ccle=dat2[,x],
                                                ccle.name=colnames(dat2)[x],
                                                pvalue))
        corr <- do.call(rbind,corr)
        out[['pearson']][[s]] <- as.data.frame(corr)
    }
    if (method == "all" | method == "manhattan"){
        man <- lapply(seq_len(ncol(dat2)),
                      function(x) pair_manhattan(cell=cell_cn,
                                                 ccle=dat2[,x],
                                                 ccle.name=colnames(dat2)[x],
                                                 pvalue))
        man <- do.call(rbind,man)
        out[['manhattan']][[s]] <- as.data.frame(man)
    }
    if (method == "all" | method == "euclidean"){
        eu <- lapply(seq_len(ncol(dat2)),
                     function(x) pair_euclidean(cell=cell_cn,
                                                ccle=dat2[,x],
                                                ccle.name=colnames(dat2)[x],
                                                pvalue))
        eu <- do.call(rbind,eu)
        out[['euclidean']][[s]] <- as.data.frame(eu)
    }
    if (method == "all" | method == "cosine"){
        cos <- lapply(seq_len(ncol(dat2)),
                     function(x) pair_cosine(cell=cell_cn,
                                             ccle=dat2[,x],
                                             ccle.name=colnames(dat2)[x],
                                             pvalue))
        cos <- do.call(rbind,cos)
        out[['cosine']][[s]] <- as.data.frame(cos)
    }

    pearson<-data.table::rbindlist(out[['pearson']], idcol=TRUE)
    manhattan<-data.table::rbindlist(out[['manhattan']], idcol=TRUE)
    euclidean<-data.table::rbindlist(out[['euclidean']], idcol=TRUE)
    cosine<-data.table::rbindlist(out[['cosine']], idcol=TRUE)

    #join all metrics and return a data.frame
    if (method=="all" & pvalue==TRUE){
        measures<-cbind(pearson,manhattan[,c(3,4)],euclidean[,c(3,4)],cosine[,c(3,4)])
        colnames(measures)[1]<-"fileid"
        return(measures)
    }
    if (method=="all" & pvalue==FALSE){
        measures<-cbind(pearson,manhattan[,3],euclidean[,3],cosine[,3])
        colnames(measures)[1]<-"fileid"
        return(measures)
    }
    #or return only the desired metric
    else if (method=="pearson"){
        colnames(pearson)[1]<-"fileid"
        return(pearson)
    }
    else if (method=="manhattan"){
        colnames(manhattan)[1]<-"fileid"
        return(manhattan)
    }
    else if (method=="euclidean"){
        colnames(euclidean)[1]<-"fileid"
        return(euclidean)
    }
    else if (method=="cosine"){
        colnames(cosine)[1]<-"fileid"
        return(cosine)
    }
}


#' @title Get top hits
#' @description This function get the top hit for each comparison per sample
#' @name getTopHit
#'
#' @param samples vector with name of samples
#' @param measure dataframe with similarity metrics
#' @param method similarity metric to use for ordering. Default is all
#' @return dataframe with similarity metrics for all comparisons
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' ccle_cn<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1:2])
#' exp_cell<-as.matrix(ccle_cn[,1])
#' colnames(exp_cell)<-unique(cells_segcn$sample)[1]
#' measures<-getSimilarities(dat1=exp_cell,dat2=ccle_cn,method="pearson")
#' tophits<-getTopHit(samples=unique(cells_segcn$sample)[1],measure=measures,method="pearson")

getTopHit<-function(samples, measure, method="all"){
    colnames(measure)[1]<-"cellid"
    tops<-list()
    if (method=="all" | method=="pearson"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(-m$r),]
            top<-m[1,c(seq_len(2))]
            tops[['pearson']][[i]]<-top
        }
    }
    if (method=="all" | method=="cosine"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(-m$cos_sim),]
            top<-m[1,c(seq_len(2))]
            tops[['cosine']][[i]]<-top
        }
    }
    if (method=="all" | method=="manhattan"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(m$manhattan),]
            top<-m[1,c(seq_len(2))]
            tops[['manhattan']][[i]]<-top
        }
    }
    if (method=="all" | method=="euclidean"){
        for (i in samples){
            m<-measure[measure$cellid%in%i,]
            m<-m[order(m$euclidean),]
            top<-m[1,c(seq_len(2))]
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
#' @importFrom methods is
#'
#' @param dat bin-level copy numbers
#' @return matrix with copy numbers of all samples
#' @export
#' @examples
#' posBins<-lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' exp_cell<-getCNbins(posBins=posBins, data=cells_segcn, samples=unique(cells_segcn$sample)[1])
#' colnames(exp_cell)<-unique(cells_segcn$sample)[1]
#' cnmatrix<-getInputmatrix(exp_cell)

getInputmatrix<-function(dat){
    #prepare matrix with cell lines copy-number
    if(is(dat, "QDNAseqCopyNumbers") == TRUE) {
        cn_matrix <- dat@assayData$copynumber
    } else {
        cn_matrix <- dat
    }
    return(cn_matrix)
}


#' @title Get list of samples included in the closest cluster to the input sample
#' @description This function cluster samples including the input sample or not.
#' If the input sample is not include, cosine similarity between the input sample
#' and the cluster centers is computed to obtain the list of samples included
#' in the closest cluster to the input sample
#' @name getClusterSamples
#'
#' @param matrix matrix with variables in rows and samples in columns
#' @param cell one-column matrix with data of input sample, with one row per variable
#' @param include a logical value to indicate if input sample should be used for clustering or not
#' @return vector with list of samples in the closest cluster of the input sample
#' @export
#' @examples
#' matrix <- cbind(s1=runif(10, min=0, max=1), s2=runif(10, min=0, max=1),s3=runif(10, min=0, max=1))
#' rownames(matrix)<-paste0("c", seq(1,10))
#' cell<-matrix(runif(n = 3, min = 0, max = 1))
#' colnames(cell)<-"input"
#' samples<-getClusterSamples(t(matrix), cell)

getClusterSamples<-function(matrix, cell, include=TRUE){
    set.seed(1)
    if (include==TRUE) {
        m<-cbind(matrix,cell)
        km.res<-akmeans::akmeans(t(m), d.metric=2, ths3=0.8, mode=3, verbose=FALSE)
        cluster<-km.res$cluster[names(km.res$cluster)==colnames(cell)]
        samples<-names(km.res$cluster)[km.res$cluster==cluster]

    } else {
        km.res<-akmeans::akmeans(t(matrix), d.metric=2, ths3=0.8, mode=3, verbose=FALSE)
        centers<-t(km.res$centers)
        colnames(centers)<-seq(1,ncol(centers))
        cos <- lapply(1:ncol(centers),
                      function(x) pair_cosine(cell=cell,
                                              ccle=centers[,x],
                                              ccle.name=colnames(centers)[x],
                                              pvalue=FALSE))
        cos <- as.data.frame(do.call(rbind,cos))
        cos$cos_sim <- as.numeric(cos$cos_sim)
        cos <- cos[order(-cos$cos_sim),]
        samples <- names(km.res$cluster)[km.res$cluster==cos[1,1]]
    }
    return(samples)
}

