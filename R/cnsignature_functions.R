##########################################
## FUNCTIONS FOR QUANTIFYING SIGNATURES ##
##########################################

#### Main functions ####

#' @title Extract features of copy number profiles
#' @description This function extract feactures of CNAs.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name extractCopynumberFeatures
#'
#' @param CN_data segment table with copy numbers
#' @param cores cores to use. Default is 1
#'
#' @return A list with features per sample
#' @examples
#' @export

extractCopynumberFeatures<-function(CN_data, cores = 1){
    #get chromosome lengths
    chrlen<-chrlen

    #get centromere locations
    gaps<-gaps
    centromeres<-gaps[gaps[,8]=="centromere",]

    if(cores > 1) {
        doMC::registerDoMC(cores)

        temp_list = foreach::foreach(i=1:6) %dopar% {
            if(i == 1){
                list(segsize = getSegsize(CN_data) )
            } else if (i == 2) {
                list(bp10MB = getBPnum(CN_data,chrlen) )
            } else if (i == 3) {
                list(osCN = getOscilation(CN_data,chrlen) )
            } else if (i == 4) {
                list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
            } else if (i == 5) {
                list(changepoint = getChangepointCN(CN_data) )
            } else {
                list(copynumber = getCN(CN_data) )
            }

        }
        unlist( temp_list, recursive = FALSE )
    } else {

        segsize<-getSegsize(CN_data)
        bp10MB<-getBPnum(CN_data,chrlen)
        osCN<-getOscilation(CN_data,chrlen)
        bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
        changepoint<-getChangepointCN(CN_data)
        copynumber<-getCN(CN_data)

        list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
    }

}

#' @title Quantifying signatures
#' @description This function quantify copy-number signatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name quantifySignatures
#'
#' @param sample_by_component matrix with components per sample
#' @param component_by_signature matrix with components per signature
#'
#' @return A matrix with signature exposures per sample
#' @examples
#' @export

quantifySignatures<-function(sample_by_component,component_by_signature=NULL){
    if(is.null(component_by_signature))
    {
        component_by_signature <- feat_sig_mat
    }
    signature_by_sample<-YAPSA::LCD(t(sample_by_component),
                                    YAPSA:::normalize_df_per_dim(component_by_signature,2))
    signature_by_sample<-normaliseMatrix(signature_by_sample)
    signature_by_sample
}

#' @title Get components per sample
#' @description This function get components per sample.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name generateSampleByComponentMatrix
#'
#' @param CN_features list with features per sample
#' @param all_components components of copy-number signatures to obtain
#' @param cores number of cores to use. Default is 1
#' @param rowIter iterations for calculation of sum of posteriors
#' @param subcores number of subcores to use. Default is 2
#'
#' @return A matrix with components per sample
#' @examples
#' @export

generateSampleByComponentMatrix<-function(CN_features, all_components=NULL, cores = 1, rowIter = 1000, subcores = 2){
    if(is.null(all_components))
    {
        all_components<-readRDS("data/component_parameters.rds")
    }

    if(cores > 1){

        feats = c( "segsize", "bp10MB", "osCN", "changepoint", "copynumber", "bpchrarm" )
        doMC::registerDoMC(cores)

        full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
            calculateSumOfPosteriors(CN_features[[feat]],all_components[[feat]],
                                     feat, rowIter = rowIter, cores = subcores)
        }
    } else {
        full_mat<-cbind(
            calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
            calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
            calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
            calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
            calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
            calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"))
    }

    rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
    full_mat[is.na(full_mat)]<-0
    full_mat
}


#### Helper functions ####

#' @title Calculate sum of posteriors
#' @description Helper function for generateSampleByComponentMatrix.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name calculateSumOfPosteriors
#'
#' @param CN_feature feature of the copy-number profiles
#' @param components components of copy-number signatures
#' @param name name of the component
#' @param rowIter iterations for calculation of sum of posteriors. Default is 1000
#' @param cores number of cores to use. Default is 1
#'
#' @return Sum of posteriors per component
#' @examples
#' @export

calculateSumOfPosteriors<-function(CN_feature,components,name, rowIter = 1000, cores = 1){
    if(cores > 1){
        len = dim(CN_feature)[1]
        iters = floor( len / rowIter )
        lastiter = iters[length(iters)]

        registerDoMC(cores)
        curr_posterior = foreach( i=0:iters, .combine=rbind) %dopar% {
            start = i*rowIter+1
            if(i != lastiter) { end = (i+1)*rowIter } else { end = len }
            flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[start:end,2])))
        }
    } else {
        curr_posterior<-flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[,2])))
    }

    mat<-cbind(CN_feature,curr_posterior)
    posterior_sum<-c()

    ## foreach and parallelising doesn't make the following code faster.
    for(i in unique(mat$ID))
    {
        posterior_sum<-rbind(posterior_sum,colSums(mat[mat$ID==i,c(-1,-2)]))
    }
    params<-flexmix::parameters(components)
    if(!is.null(nrow(params)))
    {
        posterior_sum<-posterior_sum[,order(params[1,])]
    }
    else
    {
        posterior_sum<-posterior_sum[,order(params)]
    }
    colnames(posterior_sum)<-paste0(name,1:ncol(posterior_sum))
    rownames(posterior_sum)<-rownames(unique(mat$ID))
    posterior_sum
}


#' @title Get length of segments
#' @description Helper function for extractCopynumberFeatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getSegsize
#'
#' @param abs_profiles absolute copy-number profiles
#'
#' @return Matrix with segment sizes per sample
#' @examples
#' @export

getSegsize<-function(abs_profiles){
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        seglen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
        seglen<-seglen[seglen>0]
        out<-rbind(out,cbind(ID=rep(i,length(seglen)),value=seglen))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


#' @title Get counts of brekpoints per 10Mb
#' @description Helper function for extractCopynumberFeatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getBPnum
#'
#' @param abs_profiles absolute copy-number profiles
#' @param chrlen chromosome length
#'
#' @return Matrix with counts of break points per 10Mb per sample
#' @examples
#' @export

getBPnum<-function(abs_profiles,chrlen){
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        allBPnum<-c()
        for(c in chrs)
        {
            currseg<-segTab[segTab$chromosome==c,]
            intervals<-seq(1,chrlen[chrlen[,1]==paste0("chr",c),2]+10000000,10000000)
            res <- hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
            allBPnum<-c(allBPnum,res)
        }
        out<-rbind(out,cbind(ID=rep(i,length(allBPnum)),value=allBPnum))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


#' @title Get length of chains of oscilating copy numbers
#' @description Helper function for extractCopynumberFeatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getBPnum
#'
#' @param abs_profiles absolute copy-number profiles
#' @param chrlen chromosome length
#'
#' @return Matrix with length of chains of oscilating copy numbers
#' @examples
#' @export

getOscilation<-function(abs_profiles,chrlen){
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        oscCounts<-c()
        for(c in chrs)
        {
            currseg<-segTab[segTab$chromosome==c,"segVal"]
            currseg<-round(as.numeric(currseg))
            if(length(currseg)>3)
            {
                prevval<-currseg[1]
                count=0
                for(j in 3:length(currseg))
                {
                    if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
                    {
                        count<-count+1
                    }else{
                        oscCounts<-c(oscCounts,count)
                        count=0
                    }
                    prevval<-currseg[j-1]
                }
            }
        }
        out<-rbind(out,cbind(ID=rep(i,length(oscCounts)),value=oscCounts))
        if(length(oscCounts)==0)
        {
            out<-rbind(out,cbind(ID=i,value=0))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


#' @title Get the breakpoint count per chromosome arm
#' @description Helper function for extractCopynumberFeatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getCentromereDistCounts
#'
#' @param abs_profiles absolute copy-number profiles
#' @param centromeres centromere positions
#' @param chrlen chromosome length
#'
#' @return Matrix with breakpoint count per chromosome arm
#' @examples
#' @export

getCentromereDistCounts<-function(abs_profiles,centromeres,chrlen){
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        chrs<-unique(segTab$chromosome)
        all_dists<-c()
        for(c in chrs)
        {
            if(nrow(segTab)>1)
            {
                starts<-as.numeric(segTab[segTab$chromosome==c,2])[-1]
                segstart<-as.numeric(segTab[segTab$chromosome==c,2])[1]
                ends<-as.numeric(segTab[segTab$chromosome==c,3])
                segend<-ends[length(ends)]
                ends<-ends[-length(ends)]
                centstart<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
                centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
                chrend<-chrlen[substr(chrlen[,1],4,5)==c,2]
                ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
                ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
                ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
                ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
                ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
                ndist<-apply(ndist,1,min)

                all_dists<-rbind(all_dists,sum(ndist>0))
                all_dists<-rbind(all_dists,sum(ndist<=0))
            }
        }
        if(nrow(all_dists)>0)
        {
            out<-rbind(out,cbind(ID=i,ct1=all_dists[,1]))
        }
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


#' @title Get difference in copy number between adjacent segments
#' @description Helper function for extractCopynumberFeatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getChangepointCN
#'
#' @param abs_profiles absolute copy-number profiles
#'
#' @return Matrix with difference in copy number between adjacent segments
#' @examples
#' @export

getChangepointCN<-function(abs_profiles){
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        chrs<-unique(segTab$chromosome)
        allcp<-c()
        for(c in chrs)
        {
            currseg<-as.numeric(segTab[segTab$chromosome==c,"segVal"])
            allcp<-c(allcp,abs(currseg[-1]-currseg[-length(currseg)]))
        }
        if(length(allcp)==0)
        {
            allcp<-0 #if there are no changepoints
        }
        out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=allcp))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


#' @title Get copy numbers
#' @description Helper function for extractCopynumberFeatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getCN
#'
#' @param abs_profiles absolute copy-number profiles
#'
#' @return Matrix with copy numbers
#' @examples
#' @export

getCN<-function(abs_profiles){
    out<-c()
    samps<-getSampNames(abs_profiles)
    for(i in samps)
    {
        if(class(abs_profiles)=="QDNAseqCopyNumbers")
        {
            segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
        }
        else
        {
            segTab<-abs_profiles[[i]]
            colnames(segTab)[4]<-"segVal"
        }
        segTab$segVal[as.numeric(segTab$segVal)<0]<-0
        cn<-as.numeric(segTab$segVal)
        out<-rbind(out,cbind(ID=rep(i,length(cn)),value=cn))
    }
    rownames(out)<-NULL
    data.frame(out,stringsAsFactors = F)
}


#' @title Get sample names
#' @description Helper function for getting the name of samples.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getSampNames
#'
#' @param abs_profiles absolute copy-number profiles
#'
#' @return name of each sample
#' @examples
#' @export

getSampNames<-function(abs_profiles){
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
        samps<-colnames(abs_profiles)
    }
    else
    {
        samps<-names(abs_profiles)
    }
    samps
}


#' @title Get segmented table of each sample
#' @description Helper function for getting segment table with copy numbers.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name getSegTable
#'
#' @param x each object included
#'
#' @return segment table per sample
#' @examples
#' @export

getSegTable<-function(x){
    dat<-x
    sn<-Biobase::assayDataElement(dat,"segmented")
    fd <- Biobase::fData(dat)
    fd$use -> use
    fdfiltfull<-fd[use,]
    sn<-sn[use,]
    segTable<-c()
    for(c in unique(fdfiltfull$chromosome))
    {
        snfilt<-sn[fdfiltfull$chromosome==c]
        fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
        sn.rle<-rle(snfilt)
        starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
        ends <- cumsum(sn.rle$lengths)
        lapply(1:length(sn.rle$lengths), function(s) {
            from <- fdfilt$start[starts[s]]
            to <- fdfilt$end[ends[s]]
            segValue <- sn.rle$value[s]
            c(fdfilt$chromosome[starts[s]], from, to, segValue)
        }) -> segtmp
        segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
        segTable<-rbind(segTable,segTableRaw)
    }
    colnames(segTable) <- c("chromosome", "start", "end", "segVal")
    segTable
}


#' @title Normalize matrix
#' @description Helper function for quantifySignatures.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name normaliseMatrix
#'
#' @param signature_by_sample matrix with signature levels per sample
#' @param sig_thresh threshold for each signature. Default is 0.01
#'
#' @return matrix with normalized activity of each signature per sample
#' @examples
#' @export

normaliseMatrix<-function(signature_by_sample,sig_thresh=0.01){
    norm_const<-colSums(signature_by_sample)
    sample_by_signature<-apply(signature_by_sample,1,function(x){x/norm_const})
    sample_by_signature<-apply(sample_by_signature,1,lower_norm,sig_thresh)
    signature_by_sample<-t(sample_by_signature)
    norm_const<-apply(signature_by_sample,1,sum)
    sample_by_signature<-apply(signature_by_sample,2,function(x){x/norm_const})
    signature_by_sample<-t(sample_by_signature)
    signature_by_sample
}


#' @title Get lower normalized value
#' @description Helper function for normaliseMatrix.
#' Approach developed by Geoff Macintyre et al. 2018
#' @name lower_norm
#'
#' @param x value
#' @param sig_thresh threshold for each signature. Default is 0.01
#'
#' @return lower normalized value
#' @examples
#' @export

lower_norm<-function(x,sig_thresh=0.01){
    new_x<-x
    for(i in 1:length(x))
    {
        if(x[i]<sig_thresh)
        {
            new_x[i]<-0
        }
    }
    new_x
}
