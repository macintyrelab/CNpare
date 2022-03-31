###########################################
## FUNCTIONS FOR QUANTIFYING DIFFERENCES ##
###########################################


#' @title Alignment of genomic positions of segments
#' @description This function unifies boundaries of copy-number segments.
#' Before to apply this function, please double-check that both profiles include
#' the same chromosomes
#' @name unifySegments
#'
#' @param posSeg dataframe with segments of the reference sample
#' @param data dataframe with segments of the experimental model
#' @param window boundary margin for smoothing segments. Default is 100000
#'
#' @return dataframe with segments of the model aligned to those in reference
#' @export
#' @examples
#' unify<-unifySegments(posSeg=cells_segcn[cells_segcn$sample=="OVKATE",],
#'     data=cells_segcn[cells_segcn$sample=="OV-90",])

unifySegments<-function(posSeg,data,window=100000){
    cn_out<-c()
    for(c in unique(posSeg[,1])){
        ref<-posSeg[posSeg$chromosome==c,]
        exp<-data[data$chromosome==c,]

        #unify & collapse segments within window distance
        new_seg<-getNewSegments(ref,exp,window)

        #generate output table with unified copy number
        cn<-c()
        ref_cn<-getUnifiedCN(new_seg,ref)
        exp_cn<-getUnifiedCN(new_seg,exp)
        cn<-as.data.frame(cbind(chromosome=c,start=new_seg[,1],end=new_seg[,2],segVal=ref_cn,sample=unique(data$sample),segVal_B=exp_cn))
        cn_out<-rbind(cn_out,cn)
    }
    cn_out[,c(2:4,6)]<-sapply(cn_out[,c(2:4,6)],as.numeric)
    return(cn_out)
}


#' @title Create new segments for unifying profiles in a specific chromosome
#' @description This function creates unified and collapse segments within a window distance.
#' This function in used in unifySegments function, where a table with new segments per
#' chromosome is generated for generating a unified copy number table.
#' @name getNewSegments
#'
#' @param ref dataframe with segments of the reference sample
#' @param exp dataframe with segments of the experimental model
#' @param window boundary margin for smoothing segments. Default is 100000
#'
#' @return dataframe with new segments in a specific chromosome
#' @export
#' @examples
#' posSeg=cells_segcn[cells_segcn$sample=="OVKATE",]
#' data=cells_segcn[cells_segcn$sample=="OV-90",]
#' new_seg<-getNewSegments(ref=posSeg[posSeg$chromosome==1,],
#'     exp=data[data$chromosome==1,],window=100000)

getNewSegments<-function(ref,exp,window){
    #unify chromosome boundaries in both profiles
    start<-max(min(ref[,2]),min(exp[,2]))
    end<-min(max(ref[,3]),max(exp[,3]))
    ref[1,2]<-start
    exp[1,2]<-start
    ref[nrow(ref),3]<-end
    exp[nrow(exp),3]<-end

    #unify & collapse internal segments based on reference profile
    starts<-unique(c(ref[,2],exp[,2]))
    starts<-starts[order(starts)]
    ends<-c(starts[-1]-1,end)

    #collapse segments that are within window distance
    starts_collapsed<-c()
    ends_collapsed<-c()
    if(length(starts)>1){
        for(i in 1:(length(starts)-1)){
            if(abs(starts[i+1]-starts[i])>(window+1)){
                starts_collapsed<-c(starts_collapsed,starts[i])
                ends_collapsed<-c(ends_collapsed,(starts[i+1]-1))
            }
        }
        starts_collapsed<-c(starts_collapsed,ends_collapsed[length(ends_collapsed)]+1)
        ends_collapsed<-c(ends_collapsed,ends[length(ends)])
    }else{
        starts_collapsed<-starts
        ends_collapsed<-ends
    }
    new_seg<-cbind(starts_collapsed,ends_collapsed)
    return(new_seg)
}


#' @title Get copy number values of each new unified segment from a chromosome
#' @description This function generates a vector with copy number values
#' of each new segment, which has been defined using getNewSegments function.
#' This function in used in unifySegments function, where unified copy number values
#' in each new chromosome segments is obtained to generate a unified copy number table
#' @name getUnifiedCN
#'
#' @param new_seg dataframe defining the new unified segments in a specific chromosome
#' @param data dataframe with original segments of the sample in a specific chromosome
#'
#' @return dataframe with unified segment tables
#' @export
#' @examples
#' posSeg=cells_segcn[cells_segcn$sample=="OVKATE",]
#' data=cells_segcn[cells_segcn$sample=="OV-90",]
#' new_seg<-getNewSegments(ref=posSeg[posSeg$chromosome==1,],
#'     exp=data[data$chromosome==1,],window=250000)
#' new_cn<-getUnifiedCN(new_seg,posSeg[posSeg$chromosome==1,])

getUnifiedCN<-function(new_seg,data){
    cn<-c()
    if(nrow(new_seg)>1){
        for(i in 1:(nrow(new_seg)-1)){
            start<-new_seg[i,1]
            end<-new_seg[i,2]
            ind<-intersect(as.numeric(which(data$end>=start)), as.numeric(which(data$start<=end)))
            if(length(ind)!=0){cn<-c(cn,data[ind,4])}
            if(length(ind)==0){cn<-c(cn,data[nrow(data),4])}
        }
        #add last segment
        cn<-c(cn,data[nrow(data),4])
    }else{
        cn<-c(cn,data[nrow(data),4])
    }
    return(cn)
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

    #round absolute copy-number
    unify$segVal<-round(unify$segVal)
    unify$segVal_B<-round(unify$segVal_B)

    diff<-c()
    #Identify similar and different segments
    sameSegs<-unify[which(unify$segVal==unify$segVal_B),]
    difSegs<-unify[which(unify$segVal!=unify$segVal_B),]

    #Calculate the % of difference
    sum_diff<-sum(as.numeric(difSegs$end-difSegs$start))
    diff<-(sum_diff/genome_size)*100
}
