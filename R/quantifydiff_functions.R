###########################################
## FUNCTIONS FOR QUANTIFYING DIFFERENCES ##
###########################################


#' @title Alignment of genomic positions of segments
#' @description This function unifies boundaries of copy-number segments
#' @name unifySegments
#' @importFrom data.table as.data.table setkey foverlaps
#'
#' @param posSeg dataframe with segments of the reference sample
#' @param data dataframe with segments of the experimental model
#'
#' @return dataframe with segments of the model aligned to those in reference
#' @export
#' @examples
#' unify<-unifySegments(posSeg=cells_segcn[cells_segcn$sample=="OVKATE",],
#'     data=cells_segcn[cells_segcn$sample=="OV-90",])

unifySegments<-function(posSeg,data){
    pos<-data.table::as.data.table(posSeg[,c(5,4,1:3)])
    data<-data.table::as.data.table(data[,c(5,4,1:3)])
    data.table::setkey(pos, chromosome, start, end)
    overlap <- data.table::foverlaps(data,pos)
    cn_new<-c()
    for(i in seq_len(nrow(posSeg)))
    {
        chr<-posSeg[i,1]
        start<-posSeg[i,2]
        end<-posSeg[i,3]
        curr_cn<-overlap$i.segVal
        ind<-which(overlap$chromosome==chr & overlap$start==start & overlap$end==end)
        cn_new<-c(cn_new,median(curr_cn[c(ind[1]:ind[length(ind)])],na.rm=TRUE))
    }
    cn_out<-cbind(posSeg,cn_new)
    cn_out<-data.frame(cn_out,stringsAsFactors = FALSE)
    colnames(cn_out)<-c("chromosome","start","end","segVal","sample","segVal_B")

    #round absolute copy-number
    cn_out$segVal<-round(cn_out$segVal)
    cn_out$segVal_B<-round(cn_out$segVal_B)
    #return
    cn_out
}

#' @title Calculation of ploidy status
#' @description This function calculate the overall ploidy of samples by computing
#' the weighted mean copy number values of all segments
#' @name getPloidy
#' @importFrom stats weighted.mean
#'
#' @param events dataframe with segments of the sample
#'
#' @return ploidy status of the sample
#' @export
#' @examples
#' exp_cell=cells_segcn[cells_segcn$sample=="22RV1",]
#' ploidy<-getPloidy(exp_cell)

getPloidy<-function(events){
    events$length<-events$end-events$start
    ploidy <- stats::weighted.mean(events$segVal, events$length)
    ploidy <-round(as.numeric(ploidy))
    return(ploidy)
}


#' @title Calculation of the extent of differences between two profiles
#' @description This function computes the % genome difference
#' @name getDifference
#'
#' @param events segment table with absolute copy numbers of one sample
#' @param events_2 segment table with absolute copy numbers of the other sample
#' @param method method used for calculating difference. Options are "normalized"
#' or "non-normalized" by ploidy status. Default is "non-normalized"
#'
#' @return % genome difference between two copy-number profiles
#' @export
#' @examples
#' diff<-getDifference(events=cells_segcn[cells_segcn$sample=="OVKATE",],
#'     events_2=cells_segcn[cells_segcn$sample=="OV-90",])

getDifference<- function(events, events_2, method="non-normalized"){
    #unify segments
    unify<-unifySegments(events, events_2)

    #normalization by ploidy status
    if (method == "normalized"){
        ploidy<-getPloidy(events)
        ploidy2<-getPloidy(events_2)
        ind<-ploidy>ploidy2
        if (ind == TRUE){
            unify$segVal_B<-unify$segVal_B+(max(ploidy,ploidy2)-min(ploidy,ploidy2))
        } else {
            unify$segVal<-unify$segVal+(max(ploidy,ploidy2)-min(ploidy,ploidy2))
        }
    }

    #get genome size
    genome_size<-chr_sizes
    genome_size<-sum(genome_size$length)

    diff<-c()
    #Identify similar and different segments
    sameSegs<-unify[which(unify$segVal==unify$segVal_B),]
    difSegs<-unify[which(unify$segVal!=unify$segVal_B),]

    #Calculate the % of difference
    sum_diff<-sum(as.numeric(difSegs$end-difSegs$start))
    diff<-(sum_diff/genome_size)*100
}
