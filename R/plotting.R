######################################
## FUNCTIONS FOR DATA VISUALIZATION ##
######################################


#' @title Density plot showing the distribution of genome differences
#' @description This function draws the extent of genome differences across all comparisons
#' @name plot_diffdensity
#' @import ggplot2
#'
#' @param diff dataframe with cell names 'cellid' and genome differences 'percDiff'
#'
#' @return A density plot with the distribution. The min and max value is also reported
#' @examples
#' @export

plot_diffdensity<-function(diff){
    colnames(diff)[2]<-"percDiff"
    ggplot(data = diff, aes(x=as.numeric(percDiff)))+
        geom_density(fill='firebrick4')+
        xlim(0,100)+
        labs(title= "Distribution of the genomic differences",
             subtitle = (paste0("% from ",round(min(as.numeric(diff$percDiff)),2),
                                " to ",round(max(as.numeric(diff$percDiff)),2))),
             x="% of difference",
             y="Frequency")+
        theme_minimal()
}


#' @title Visualization of two copy-number profiles
#' @description This function draws two copy-number profiles for visual comparison
#' @name CNPlot_events
#' @import graphics
#'
#' @param events segment table with absolute copy numbers of one sample
#' @param events_2 segment table with absolute copy numbers of the other sample
#' @param plot_diff logical indicating if only segments that are different in
#' the second profile may be plotted. Default is FALSE
#' @return A plot with copy-number profiles of two samples. The % genome difference
#' between profiles is also reported.
#' @examples
#' @export

CNPlot_events <- function(events,events_2,plot_diff=FALSE){
    #calculate % differences
    unify<-unifySegments(events,events_2)
    percentage_diff<-getDifference(unify)

    #data for plotting
    chr_sizes <- read.csv(file="data/hg19.chrom.sizes.txt",sep="\t",stringsAsFactors = F,
                          header=FALSE,col.names = c("chr","length"))[1:23,]
    chr_sizes$length <- 1e-6 * chr_sizes$length #Convert sizes to Mb
    chr_sizes$offset <- cumsum(chr_sizes$length) - chr_sizes$length #Calculate the interval range of each chr

    #plotting
    if (plot_diff==FALSE){
        events<-CNconvert(events,chr_sizes)
        events_2<-CNconvert(events_2,chr_sizes)

        #plot
        plot.new()
        plot.window(xlim=c(0,sum(events$length)),ylim=c(0,10),yaxs="i",xaxs="r")
        abline(v=chr_sizes$offset,lty=2,lwd=1)
        abline(v=0,lty=1,lwd=1)
        #title
        title(main="Copy Number Profile",xlab="Chromosome",ylab="Copy Number",
              cex.lab=1.1,cex.main=1.5,font.lab=2,family="Palatino",outer=F,
              sub=paste0(events$sample," profile is ",round(percentage_diff,2)," % different to ",events_2$sample," profile"))
    } else {
        #get unified profiles
        s<-unique(events_2$sample)
        events<-unify[,c(1:5)]
        events_2<-unify[,c(1:3,6,5)]
        colnames(events_2)[4]<-"segVal"
        events_2<-events_2[which(events$segVal!=events_2$segVal),]
        events<-CNconvert(events,chr_sizes)
        events_2<-CNconvert(events_2,chr_sizes)

        #plot
        plot.new()
        plot.window(xlim=c(0,sum(events$length)),ylim=c(0,10),yaxs="i",xaxs="r")
        abline(v=chr_sizes$offset,lty=2,lwd=1)
        abline(v=0,lty=1,lwd=1)
        #title
        title(main="Copy Number Profile",xlab="Chromosome",ylab="Copy Number",
              cex.lab=1.1,cex.main=1.5,font.lab=2,family="Palatino",outer=F,
              sub=paste0(events$sample," profile is ",round(percentage_diff,2)," % different to ",s," profile"))
        #legend
        legend(x = "topright", legend = c(paste0(unique(events$sample)," profile"),paste0("Different events in ",s)),
               fill = c("red", "royalblue"), cex=0.5, xpd = TRUE, inset = c(-0.05,-0.1)) #inset = c(-0.4, -0.2)
    }

    # Set the axis
    axis(1,labels = c(1:22,"X"),tick=T,
         at = chr_sizes$offset+chr_sizes$length/2,cex.axis=1)

    lablist<-as.vector(c(0:10))
    text(par("usr")[4], seq(0, 10, by=1), labels = lablist, srt = 0, pos = 2, xpd = TRUE, offset = 1 ,
         cex=1.3)

    # Plot the first data set
    for (i in 1:nrow(events)){
        e <- events[i,]
        segments(e$start,e$segVal,e$end,e$segVal,lwd=5.,col="red")
    }

    # Plot the second data set
    for (i in 1:nrow(events_2)){
        e <- events_2[i,]
        segments(e$start,e$segVal,e$end,e$segVal,lwd=5.,col="royalblue")
    }
}


#' @title A helper function for CNPlot_events
#' @description This function prepare profiles for plotting
#' @name CNconvert
#'
#' @param e segment table with absolute copy numbers
#' @param base sizes of chromosomes
#'
#' @return segment table prepared for plotting
#' @examples
#' @export

CNconvert <- function(e,base){
    e$start <- as.numeric(e$start)
    e$end <- as.numeric(e$end)
    e$segVal <- as.numeric(e$segVal)
    e$start <- 1e-6 * e$start
    e$end <- 1e-6 * e$end
    for (c in unique(e$chromosome)){
        e$start[e$chromosome == c] <- e$start[e$chromosome == c] +
            base$offset[base$chr == c]
        e$end[e$chromosome == c]   <- e$end[e$chromosome == c]  +
            base$offset[base$chr == c]
    }
    e$length <- e$end - e$start
    return(e)
}


#' @title Plots a PCA reduction colored by cluster group
#' @description This function draws a cluster plot with samples grouped by
#' exposure levels to ovarian copy-number signatures
#' @name plotClusters
#' @import factoextra
#'
#' @param signs matrix with levels of exposure to each signature
#' @param palette vector with colors to use for each group. The length of this
#' vector must be equal to the number of cluster (k)
#' @param k number of clusters. It is recommended to estimate the
#' optimal number of clusters before
#' @return A plot with samples clustered by similarity in signature exposition
#' @examples
#' @export

plotClusters<-function(signs, palette, k){
    if (length(palette) != k) {
        stop('Number of colors must be equal to number of clusters')
    }
    km.res <- kmeans(signs, k, nstart=25)
    p<-fviz_cluster(km.res, signs,
                    palette = palette,
                    ellipse.type = "euclid", # Concentration ellipse
                    star.plot = TRUE, # Add segments from centroids to items
                    repel = TRUE, # Avoid label overplotting (slow)
                    ggtheme = theme_minimal())
    p
}
