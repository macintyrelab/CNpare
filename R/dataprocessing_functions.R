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
#' cell_bin <- getCNbins(posBins=posBins,
#'     data=cells_segcn[cells_segcn$sample=="22RV1",],
#'     samples="22RV1")

getCNbins <- function(posBins,data,samples){
    pb=data.table::rbindlist(posBins)
    CNmatrix <- matrix(nrow = nrow(pb), ncol = length(samples)) #Create the matrix
    colnames(CNmatrix) <- samples

    for (b in seq_len(nrow(pb))){
        chrom <- as.character(pb[b,1])
        start <- as.numeric(pb[b,2])
        end   <- as.numeric(pb[b,3])
        cn <- data[(data$chromosome %in% chrom & data$start<=start & data$end>=end), ]

        for (s in seq_len(length(samples))){
            if (nrow(cn)!=0){
                segVal <- cn[cn$sample==samples[s], "segVal"]
                CNmatrix[b,s] <- ifelse(length(segVal)!=0, segVal, NA)
            } else {
                CNmatrix[b,s] <- NA
            }
        }
    }
    return(CNmatrix)
}

#### Helper functions ####

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
