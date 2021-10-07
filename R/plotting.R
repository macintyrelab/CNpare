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
#' @export
#' @examples
#' differences <- as.data.frame(cbind(cellid=c(1:5), diff=runif(5, min=0, max=100)))
#' plot_diffdensity(differences)

plot_diffdensity<-function(diff){
    colnames(diff)[2]<-"percDiff"
    ggplot(data=diff, aes(x=as.numeric(percDiff)))+
        geom_density(fill='firebrick4')+
        xlim(0,100)+
        labs(title="Distribution of the genomic differences",
            subtitle=(paste0("% from ", round(min(as.numeric(diff$percDiff)), 2),
                                " to ", round(max(as.numeric(diff$percDiff)), 2))),
            x="Extent of genome difference (%)",
            y="Frequency")+
        theme_minimal()
}


#' @title Visualization of two copy-number profiles
#' @description This function draws two copy-number profiles for visual comparison
#' @name CNPlot_events
#' @importFrom graphics abline axis hist legend par plot.new plot.window segments text title
#'
#' @param events segment table with absolute copy numbers of one sample
#' @param events_2 segment table with absolute copy numbers of the other sample
#' @param method_diff method used for calculating difference. Options are "normalized" or
#' "non-normalized" by ploidy status
#' @param plot_diff logical indicating if only segments that are different in
#' the second profile may be plotted. Default is FALSE
#' @return A plot with copy-number profiles of two samples. The % genome difference
#' between profiles is also reported.
#' @export
#' @examples
#' exp_cell=cells_segcn[cells_segcn$sample=="22RV1",]
#' mod_cell=cells_segcn[cells_segcn$sample=="A172",]
#' CNPlot_events(exp_cell,mod_cell, method_diff="non-normalized")

CNPlot_events <- function(events, events_2, method_diff, plot_diff=FALSE){
    ##Calculate % differences
    percentage_diff<-getDifference(events, events_2, method_diff)

    ##Plotting
    #data for plotting
    chr_sizes <- chr_sizes
    chr_sizes$length <- 1e-6 * chr_sizes$length #Convert sizes to Mb
    chr_sizes$offset <- cumsum(chr_sizes$length) - chr_sizes$length #Calculate the interval range of each chr

    if (plot_diff == FALSE){
        events<-CNconvert(events,chr_sizes)
        events_2<-CNconvert(events_2,chr_sizes)

        #plot
        plot.new()
        plot.window(xlim=c(0,sum(events$length)),ylim=c(0,10),yaxs="i",xaxs="r")
        abline(v=chr_sizes$offset,lty=2,lwd=1)
        abline(v=0,lty=1,lwd=1)

        #title
        title(main="Copy Number Profile",xlab="Chromosome",ylab="Copy Number",
            cex.lab=1.1,cex.main=1.5,font.lab=2,outer=FALSE,
            sub=paste0(events$sample," profile is ",round(percentage_diff,2)," % different to ",events_2$sample," profile"))
    } else {
        #get unified profiles
        unify<-unifySegments(events, events_2)
        s<-unique(events_2$sample)
        events<-unify[,c(seq_len(5))]
        events_2<-unify[,c(seq_len(3),6,5)]
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
            cex.lab=1.1,cex.main=1.5, font.lab=2, outer=FALSE,
            sub=paste0(events$sample," profile is ",round(percentage_diff,2)," % different to ",s," profile"))
        #legend
        legend(x = "topright", legend = c(paste0(unique(events$sample)," profile"),paste0("Different events in ",s)),
            fill = c("red", "royalblue"), cex=0.5, xpd=TRUE, inset=c(-0.05,-0.1)) #inset = c(-0.4, -0.2)
    }

    # Set the axis
    axis(1, labels=c(seq_len(22),"X"), tick=TRUE,
        at = chr_sizes$offset+chr_sizes$length/2,cex.axis=1)

    lablist<-as.vector(c(0:10))
    text(par("usr")[4], seq(0, 10, by=1), labels=lablist, srt=0, pos=2,
        xpd=TRUE, offset=1, cex=1.3)

    # Plot the first data set
    for (i in seq_len(nrow(events))){
        e <- events[i,]
        segments(e$start,e$segVal,e$end,e$segVal,lwd=5.,col="red")
    }

    # Plot the second data set
    for (i in seq_len(nrow(events_2))){
        e <- events_2[i,]
        segments(e$start,e$segVal,e$end,e$segVal,lwd=5.,col="royalblue")
    }
}


#' @title A helper function for CNPlot_events
#' @description This function prepare profiles for plotting
#'
#' @param e segment table with absolute copy numbers
#' @param base sizes of chromosomes
#'
#' @return segment table prepared for plotting
#' @export
#' @examples
#' exp_cell=cells_segcn[cells_segcn$sample=="22RV1",]
#' base=CNpare:::chr_sizes
#' base$length <- 1e-6 * base$length #Convert sizes to Mb
#' base$offset <- cumsum(base$length) - base$length #Calculate the interval range of each chr
#' events<-CNconvert(exp_cell,base)

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


#' @title Plots a scatter plots with samples clustered
#' @description This function draws a scatter plot with samples grouped based of signature
#' exposures. Only the closest cluster to the input sample is colored.
#' @name plotClusters
#' @importFrom stats var
#'
#' @param matrix matrix variables in rows and samples in columns
#' @param samples vector with samples in the closest cluster of the input sample
#' @return A plot with samples clustered by similarity in signature exposition
#' @export
#' @examples
#' matrix <- cbind(s1=runif(10, min=0, max=1), s2=runif(10, min=0, max=1),
#'     s3=runif(10, min=0, max=1))
#' rownames(matrix)<-paste0("c", seq(1,10))
#' cell<-matrix(runif(n = 3, min = 0, max = 1))
#' samples<-getClusterSamples(t(matrix), cell, include=FALSE)
#' plotClusters(t(matrix), samples)


plotClusters<-function(matrix, samples){
    #perform clustering
    set.seed(1)
    km.res<-akmeans::akmeans(t(matrix), d.metric=2, ths3=0.8, mode=3, verbose=FALSE)
    #cluster where our sample is included (sample list from getClusterSamples)
    cluster<-unique(km.res$cluster[names(km.res$cluster) %in% samples])

    #get center values of clusters and select signatures with highest variability among centers
    centers<-as.data.frame(t(km.res$centers))
    centers$var<-apply(centers, 1, var)
    centers<-centers[order(-centers$var),]
    sig<-paste0("s",rownames(centers[c(1:2),]))

    #plot
    plot.data<-as.data.frame(t(matrix[rownames(matrix)%in%sig,]))
    x=colnames(plot.data)[1]
    y=colnames(plot.data)[2]
    plot.data$cluster<-as.factor(km.res$cluster)
    plot.data$sample<-rownames(plot.data)

    col<-matrix(nrow=length(unique(plot.data$cluster)), ncol=1)
    col[,1]<-"gray70"
    col[cluster,1]<-"#00AFBB"

    p<-ggpubr::ggscatter(plot.data, x, y,
                         color="cluster",
                         shape="cluster",
                         palette=col[,1],
                         ellipse = TRUE, ellipse.alpha =  0.2,
                         ellipse.type = "convex", ellipse.level = 0.95,
                         mean.point = TRUE,
                         size = 1.5,
                         star.plot = TRUE,
                         label="sample", repel=TRUE,
                         ggtheme=ggplot2::theme_minimal())
    p
}
