###################################
## FUNCTIONS FOR DATA PROCESSING ##
###################################

#### Main functions ####

#' @title Get genomic positions of bins
#' @description This function split genome into evenly sized bins
#' @name getBinsStartsEnds
#'
#' @param window numeric variable with the bin size. Default is 500kb
#' @param chr character variable with the chromosome name
#' @param lengthChr numeric variable with the chromosome length
#'
#' @return list with genomic positions of each bin
#' @export
#' @examples
#' posBins <- lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))

getBinsStartsEnds <- function(window=500000,chr,lengthChr){
    divideChr <- seq(0, lengthChr, window)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]

    binsSe <- list(chromosome=chr,starts=starts,ends=ends)
    binsSe
}


#' @title Get profiles with detectable CIN
#' @description This function select samples with detectable CIN
#' @name getCINProfiles
#'
#' @param segcn segment table of copy numbers of all samples
#' @param samples vector with sample names
#'
#' @return segment table of copy numbers of samples with detectable CIN
#' @export
#' @examples
#' cin.profiles <- getCINProfiles(segcn=cells_segcn,
#'     samples=unique(cells_segcn$sample)[1:10])

getCINProfiles <- function(segcn,samples){
    profiles <- getProfiles(segcn, samples)
    rounded_profiles <- getSegRounded(profiles)
    cin.samples <- getCINSamp(profiles=rounded_profiles)
    cin_segcn <- segcn[segcn$sample%in%cin.samples,]
    return(cin_segcn)
}


#' @title Get copy-number per bins
#' @description This function transform segment tables to bin tables
#' @name getCNbins
#'
#' @param posBins list with genomic positions of bins. Each list contains
#' data from each chromosome. Obtained from getBinsStartsEnds
#' @param data segment table of copy numbers of all samples
#' @param samples vector with sample names
#'
#' @return bin table of copy numbers of all samples
#' @export
#' @examples
#' posBins <- lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' cells_bin <- getCNbins(posBins=posBins,
#'     data=cells_segcn[cells_segcn$sample=="22RV1",], samples="22RV1")

getCNbins<- function(posBins,data,samples){
    pb=data.table::rbindlist(posBins)
    out<-list()
    out<-lapply(seq_len(nrow(pb)), function(b) getCNbins.bin(b,pb,data,samples))
    cn<-as.matrix(do.call(rbind,out))
    colnames(cn)<-samples
    return(cn)
}

#### Helper functions ####

#' @title Get copy-number per bins in a sample
#' @description This is a helper function for transforming segment tables to bin tables
#' @name getCNbins.sample
#'
#' @param s sample name
#' @param cn matrix with copy numbers in a bin (row) in samples (columns)
#'
#' @return a vector with copy numbers of samples in one bin
#' @export
#' @examples
#' posBins <- lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' pb=data.table::rbindlist(posBins)[1:20]
#' samp=unique(cells_segcn$sample)[1:2]
#' data=cells_segcn[cells_segcn$sample%in%samp,]
#' chrom<-as.character(pb[2,1])
#' start<-as.numeric(pb[2,2])
#' end<-as.numeric(pb[2,3])
#' cn<-data[(data$chromosome%in%chrom & data$start<=start & data$end>=end),]
#' cell_bin <- getCNbins.sample(s=samp[1],cn)

getCNbins.sample <- function(s,cn){
    if (nrow(cn)!=0){
        segVal <- cn[cn$sample==s, "segVal"]
        segVal <- ifelse(length(segVal)!=0, segVal, NA)
    } else {
        segVal <- NA
    }
    return(segVal)
}

#' @title Get copy-number in a bin
#' @description This is a helper function for transforming segment tables to bin tables
#' @name getCNbins.bin
#'
#' @param b bin number
#' @param pb matrix with genomic positions of bins
#' @param data segment table of copy numbers of one bin
#' @param samples vector with sample names
#'
#' @return a list of copy-number values per bin
#' @export
#' @examples
#' posBins <- lapply(seq_len(22),function(chr)
#'     getBinsStartsEnds(window=500000, chr, lengthChr[chr]))
#' samples=unique(cells_segcn$sample)[1:2]
#' data=cells_segcn[cells_segcn$sample%in%samples,]
#' cn_bin <- getCNbins.bin(b=2,pb=data.table::rbindlist(posBins)[1:20],data,samples)

getCNbins.bin <- function(b,pb,data,samples){
    out<-list()
    chrom<-as.character(pb[b,1])
    start<-as.numeric(pb[b,2])
    end<-as.numeric(pb[b,3])
    cn<-data[(data$chromosome%in%chrom & data$start<=start & data$end>=end),]
    out<-lapply(samples, function(s) getCNbins.sample(s,cn))
    out<-do.call(cbind,out)
    colnames(out)<-samples
    return(out)
}


#' @title Add profiles in a list
#' @description This function includes all profiles in a list
#' @name getProfiles
#'
#' @param segcn segment table of copy numbers of all samples
#' @param samples vector with sample names
#'
#' @return list of profiles
#' @export
#' @examples
#' profiles <- getProfiles(segcn=cells_segcn,
#'     samples=unique(cells_segcn$sample))

getProfiles<-function(segcn,samples){
    profiles <- list()
    for (i in samples){
        profile <- segcn[segcn$sample==i,]

        profile  <- profile[,c(seq_len(4))]
        profiles[[i]] <- profile
    }
    names(profiles) <- samples
    return(profiles)
}


#' @title Round copy-number values of segments
#' @description This function rounds copy numbers to the closest integer
#' @name getSegRounded
#'
#' @param profiles list of profiles
#'
#' @return list of profiles with copy numbers rounded
#' @export
#' @examples
#' profiles <- getProfiles(segcn=cells_segcn,
#'     samples=unique(cells_segcn$sample)[1:10])
#' rounded.profiles<-getSegRounded(profiles)

getSegRounded <- function(profiles){
    rounded_profiles<-list()
    for (i in names(profiles)){
        profile <- profiles[[i]]
        profile$segVal <- round(profile$segVal)

        segTable<-c()
        for(c in unique(profile$chromosome))
        {
            t<-profile[profile$chromosome==c,]
            rle<-rle(t$segVal)
            starts <- cumsum(c(1,rle$lengths[-length(rle$lengths)]))
            ends <- cumsum(rle$lengths)
            lapply(seq_len(length(rle$lengths)), function(s) {
                from <- t$start[starts[s]]
                to <- t$end[ends[s]]
                segValue <- rle$value[s]
                c(t$chromosome[starts[s]], from, to, segValue)
            }) -> segtmp
            segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=TRUE),stringsAsFactors=FALSE)
            segTable<-rbind(segTable,segTableRaw)
        }
        colnames(segTable) <- c("chromosome", "start", "end", "segVal")
        rounded_profiles[[i]] <- segTable
    }
    return(rounded_profiles)
}


#' @title elect samples with CIN
#' @description This function filters out samples with <20 CNAs
#' @name getCINSamp
#'
#' @param profiles list of profiles
#'
#' @return vector with sample names
#' @export
#' @examples
#' profiles <- getProfiles(segcn=cells_segcn,
#'     samples=unique(cells_segcn$sample)[1:10])
#' cin.samples<-getCINSamp(profiles)

getCINSamp<-function(profiles){
    samples<-c()
    for (i in seq_len(length(profiles))){
        profile <- profiles[[i]]
        profile <- profile[profile$chromosome != "X",]
        s<-names(profiles)[i]
        profile$dip<-profile$segVal!=2
        ndip<-nrow(profile[profile$dip,])
        if(ndip>20){
            samples<-c(samples, s)
        }
    }
    return(samples)
}
