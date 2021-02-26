Assignment 5
================

-   [Assignment Overview](#assignment-overview)
    -   [0. Getting Started](#getting-started)
    -   [1. ChIP signal tracks](#chip-signal-tracks)
    -   [2. Narrow vs Broad peaks](#narrow-vs-broad-peaks)
    -   [3. Peak calling](#peak-calling)
        -   [a. Peak calling with macs2](#a.-peak-calling-with-macs2)
        -   [b. Understanding the peaks](#b.-understanding-the-peaks)
    -   [c. Peak enrichments](#c.-peak-enrichments)

# Assignment Overview

By now you must have become vaguely familiar with ChIP-seq data but
might be a little confused about what to do after the alignment. Well,
this assignment’s aim is to walk you through a ChIP-seq pipeline
post-alignment. We will be analyzing 3 different histone modification
marks (H3K27me3, H3K4me3 and H3K27ac). In order to identify the
enrichments for each epigenetic mark, we also need to use the *input*
which represents the DNA content of the sheared chromatin sample prior
to immunoprecipitation. All the files can be found under the following
path: **/usr/local/share/data/assignment\_5/** .

-   H3K27me3 (H3K27me3\_chr3\_subset.bam)

-   H3K4me3 (H3K4me3\_chr3\_subset.bam)

-   H3K27ac (H3K27ac\_chr3\_subset.bam)

-   input (input\_chr3\_subset.bam)

A couple of things to remember:

-   When uploading your completed assignment to your GitHub directory,
    remember to specify a **github\_document** instead of the default
    *html\_document* on your .Rmd file.

-   Double check that all the files have been uploaded to your
    repository and that you are able to visualize the html-like version
    on GitHub.

-   Be mindful of your space on the server! Delete ALL unnecessary files
    and keep a clean environment.

## 0. Getting Started

We will be using a couple of new tools this time. Before we move on to
the practical part, make sure you have them all installed.

-   Integrative Genomics Viewer (IGV): Interactive tool to visualize
    different data types of genetic information (e.g. bam, bed files).
    You will install this tool to your **local computer**. To visualize
    where the reads of our ChIP analysis mapped in the genome. To
    install it, follow the instructions on this website:
    [*https://software.broadinstitute.org/software/igv/home*](https://software.broadinstitute.org/software/igv/home)

-   Deeptools
    (<https://deeptools.readthedocs.io/en/develop/index.html>): Software
    to analyze high-throughput data that allows to create easy to
    visualize figures. This will be installed on the server as part of
    your conda environment.

-   macs2
    (<https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md>):
    Tool to capture enrichments of ChIP-seq analysis. This will be
    installed on the server as part of your conda environment.

``` bash
#?# Add macs2 and deeptools to your environment created on A1 - 1 pt

#ANSWER: First, we must set up server, log into conda environment, and create a directory for this assignment. Then, we can install the tools in our environment.
ssh bkirsh@137.82.55.186
conda activate test_env
mkdir assignment5
cd assignment5
conda install -c bioconda deeptools macs2

## Install IGV to your local computer after downloading it from: https://software.broadinstitute.org/software/igv/home

#ANSWER: Done.
```

## 1. ChIP signal tracks

ChIP-seq experiments require a control against which to compare the ChIP
signnal. Often this is the “input DNA”, sheared chromatin that has not
been immmunoprecipitated. For this assignment we are using an input and
three different epigenetic marks. These histone modifications mark
states of active (H3K27ac, H3K4me3) or inactive (H3K27me3) gene
transcription and have different coverage of the genomic region where
they are located. To better visualize the differences, we will create
bigWig files from previously aligned, filtered and indexed bam files.
BigWig files are indexed, compressed files that can be used to visualize
signals across the genome. Here, we will be using it to graph the
coverage across the genome.

``` bash
## Use the bamCoverage command (included in deepTools) to convert the bam files outlined below (located in the /usr/local/share/data/assignment_5/ directory) into bigwig files that represent the read coverage (e.g. density of reads across the genome). Each file has already been indexed using sambamba and hence, has it's bam index (.bai) file associated that you will need to run bamCoverage. 
## Tip: Remember to add the .bw extension to your specified output! 

#?# Type the commands you use for converting each of the four bam files to bigwig files - 2 pts

bamCoverage -b /usr/local/share/data/assignment_5/H3K27me3_chr3_subset.bam -o H3K27me3_chr3_coverage.bw -of bigwig
bamCoverage -b /usr/local/share/data/assignment_5/H3K4me3_chr3_subset.bam -o H3K4me3_chr3_coverage.bw -of bigwig
bamCoverage -b /usr/local/share/data/assignment_5/H3K27ac_chr3_subset.bam -o H3K27ac_chr3_coverage.bw -of bigwig
bamCoverage -b /usr/local/share/data/assignment_5/input_chr3_subset.bam -o input_chr3_coverage.bw -of bigwig

### Follow this steps: 

## 1. Copy the bigwig files from the previous steps to your computer 

scp -r bkirsh@137.82.55.186:/home/bkirsh/assignment5 /mnt/c/Users/bkirs/Documents

## 2.Load all the files bigwig signal track files onto IGV on your local computer, select the "autoscale" option for each file on their individual tracks
## Tip, use the "File" tab to "load from file" option to choose the files from your computer directories

#ANSWER: Done.

## 3. Change the visualization of the files to see the following position: ** chr3:93,470,124-93,471,058 **
#?# Include a screenshot of your IGV session right after this code chunk (using the Rmarkdown syntax)  - 2 pt

#ANSWER: Done, please find the screenshot below.

#?# Explore this region by zooming in and out and navigating around. What do you see? Is there something similar across all the different files that stands out on this region? Is there anything peculiar about the DNA sequence at this locus?- 3 pt
## Tip: sometimes track signal will be truncated to a pre-set maximum. If you right-click the track label (left) and select "autoscale", it will automatically adust the scales of the track to the range of your data. Try playing around with the track settings to see their effect.

#ANSWER: Visualization with the "autoscale" option shows that the files have a similar signal distribution pattern across the region, in spite of their differing data ranges. By zooming out, we can further observe that this peak spans ~600 bps across all files, and that the nearest genes are over 500 kb away from this region. Lastly, we can see that this region is located near a centromere, where assembly has previously been shown to be difficult and result in erroneous signal. This is a characteristic location of blacklisted regions.

## 4. This file (/usr/local/share/data/assignment_5/hg38_blacklisted_regions.bed) contains the hg38 blacklisted regions. Load it into IGV along with all the other files. 

scp bkirsh@137.82.55.186:/usr/local/share/data/assignment_5/hg38.blacklist.bed /mnt/c/Users/bkirs/Documents

## 5. Look at the following region again (chr3:93,470,124-93,471,058). 
#?# What can you say now about the similarities between the files at this region? In your own words, explain what a blacklisted region is and if you think this region should be excluded a ChIP-seq analysis. - 1.5 pt

#ANSWER: Based on the overlapping nature of all five tracks, the similarity between the files at this region is likely explained by a signal that arises from a blacklisted region rather than true enrichment across all marks. In this context, blacklisted regions correspond to “hyper-chippable” regions of the genome with inexplicably enriched signal across many ChIP experiment types, regardless of what is being IPed (e.g., different types of histone marks). These false-positive peaks are atrifcats introduced by poor annotation or mapping, common at hard-masked telomeric and pericentromeric regions, rather than truly enriched regions. In my opinion, these regions should be excluded in ChIP-seq analysis because they can confound peak calling, making it hard to separate biological noise from meaningful signal. Removing such elements can help us avoid making erroneous correlations that lead to faulty biological conclusions.
```

**Add screenshot of your IGV session here:**

![Visualization of signals using IGV
(chr3:93,470,124-93,471,058).](/Users/bkirs/Documents/School/BMEG400E/GNMINFRMTCS/static/a5/igv_screenshot1.png)

## 2. Narrow vs Broad peaks

While exploring the bigwig files of the epigenetic marks on IGV, you
probably noticed that they can look very different from each other and
some of them ressemble the input more closely than others. Different
epigenetic marks can have very different signatures based on their
distribution patterns across the genome. When calling peaks, we often
lump these into two different categories of reads: broad and narrow.
Active transcription marks (H3K4me3 and H3K27ac) tend to form a sharper
coverage peaks at transcription start sites (H3K27ac is also at
enhancers), while repression marks cover a broader area (H3K27me3).
Here, we’re going to inspect their distributions relative to genes.

``` bash

## Here, we have created three bigWig track files, one for each epigenetic mark, which show the read coverage normalized using the input. They are found here: /usr/local/share/data/assignment_5/
# H3K4me3_norm.bw
# H3K27me3_norm.bw
# H3K27ac_norm.bw


## The deepTools ** computeMatrix reference-point ** command calculates scores to represent the reads mapping to specified regions of the genome across different files. 
## Use computeMatrix to compute a matrix for the signal tracks for each histone modification outlined above (which we will use to create a plot in the following step), with the following criteria: 

## - We will use the regions in reference_genes.bed located under the /usr/local/share/data/assignment_5/ directory as the basis for the plot.
## - Include the surrounding 1kb
## - Use all 3 input-normalized bigWig files (H3K4me3, H3K27ac, H3K27me3) as signal tracks
#?# Write the command you used to run it below: - 1.5 pt

computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K27me3_norm.bw /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw -R /usr/local/share/data/assignment_5/reference_genes.bed -a=1000 -b=1000 -o matrix_norm.gz

## Now that the scores matrix has been computed, we can use it to create a heatmap to provide a better visual representation of the reads distrubution across our reference genes (provided in the reference_gened.bed file)
## Use the deepTools ** plotHeatmap ** function to create a heatmap following this criteria: 
## - Use the matrix from the previous point
## - Use the Blues colormap
## - Create 3 clusters within the heatmap according to the patterns of the reads distrubution across the files using heirarchical clustering
#?# Type the command you used to run it below: - 1.5
#?# Add a screenshot of the plot right after this code chunk using Rmarkdown syntaxis - 1 pt 

plotHeatmap -m matrix_norm.gz -o matrix_norm.png --hclust 3 --colorMap Blues

#?# Explain what you are looking at (Axes, colours, curves). Where are the marks located? What are the differences between the clusters? - 3 pts

#ANSWER: First, we must copy the heatmap .png file to our local computer to visualize the file.
scp bkirsh@137.82.55.186:/home/bkirsh/assignment5/matrix_norm.png /mnt/c/Users/bkirs/Documents

#ANSWER: Using hierarchical clustering, the set of reference genes contained in the reference_genes.bed file is split into three groups based on their similarities in enrichment pattern for each of the signal files (i.e., H3K27me3, H3K4me3, and H3K27ac). The top profile plots depict the read density for each mark across the transcription start sites of all reference genes. In other words, they display the ChIP signal profile across +/- 1000 bp windows around the TSS and TES of genes for each histone modification. In these plots, the gene clusters are represented by different coloured lines (i.e., dark blue, light blue, and yellow), while the gene region and the read desnity values in the bigWig files (i.e., ChIP-seq signal intensity) correspond to the x- and y-axis, respectively. Below, the heatmaps communicate the same information, whereby the variation in colour hue reflects the magnitude of enrichment or depletion for the different marks. In this case, the darker the blue, the stronger the signal for the histone mark indicated above. Once again, the gene region is represented by the x-axis and a reference colour scale can be used to estimate signal intensity from heatmap colour values.
```

**Add screenshot here:**

![Heatmap of the distribution of input-normalized epigenetic marks
(H3K27me3, H3K4me3, and H3K27ac) relative to a set of reference
genes.](/Users/bkirs/Documents/matrix_norm.png)

``` bash
## Now the above heatmap was made with the ratio of ChIP to input. Repeat the process above, but this time using the raw bigwig files (not input-normalized). 
#?# Type the commands you used for this below - 1 pt

computeMatrix scale-regions -S H3K27me3_chr3_coverage.bw H3K4me3_chr3_coverage.bw H3K27ac_chr3_coverage.bw -R /usr/local/share/data/assignment_5/reference_genes.bed -a=1000 -b=1000 -o matrix_raw.gz

plotHeatmap -m matrix_raw.gz -o matrix_raw.png --hclust 3 --colorMap Blues

scp bkirsh@137.82.55.186:/home/bkirsh/assignment5/matrix_raw.png /mnt/c/Users/bkirs/Documents

#?# Include a screenshot of this analysis, below this code block. - 1 pt
#?# How does this compare to the input-normalized data? Why do you think this is? - 1 pt

#ANSWER: The range of the scale for the heatmaps generated with the raw ChIP-seq data is noticably larger, and the signal intensities appear to be much more binary (i.e., "all-or-nothing", white vs. dark-blue). Furthermore, the size of the clusters is different, with cluister_1 being much smaller and cluster_3 being much larger in the raw heat maps. This is because the first set of heatmaps was created by normalizing the data based on the ratio of ChIP to input values, therefore compensating for differences in sequencing depth and mapping efficiency. In this sense, the input DNA acts as a control sample that is used to define true peaks as the regions of the genome with statistically significant higher number of reads than the control/background. Since histone modifications have different shapes, levels of enrichment, and coverage of the genomic region where they are located, normalization against a background must be applied to preserve the shape of each variable’s distribution and make them easily comparable on the same scale. Otherwise, examining the global reads distribution will not take into account the differences in shape of the enrichment profile expected with different ChIP-seq targets.
```

**Add screenshot here:**

![Heatmap of the distribution of raw epigenetic marks (H3K27me3,
H3K4me3, and H3K27ac) relative to a set of reference
genes.](/Users/bkirs/Documents/matrix_raw.png)

## 3. Peak calling

Now we want to identify enriched regions of the genome for each of our
three histone marks. In order to get the enrichments, we will run the
**macs2** program to call the peaks for each epigenetic mark.

### a. Peak calling with macs2

``` bash
## Tip: Make sure to read the documentation (using the -h flag) for the *masc2 callpeak* command
## Run the callpeak command of the macs2 tool, once for each of H3K27ac, H3K27Me3, H3K4Me3
## Each time, use the input as the control 
## For H3K27Me3, you should call "broad" peaks because these signals tend to form longer domains and not acute peaks.
#?# Type the commands you used below: - 1.5 pt

macs2 callpeak -t /usr/local/share/data/assignment_5/H3K27me3_chr3_subset.bam -c /usr/local/share/data/assignment_5/input_chr3_subset.bam -f BAM -n H3K27me3 --broad
macs2 callpeak -t /usr/local/share/data/assignment_5/H3K4me3_chr3_subset.bam -c /usr/local/share/data/assignment_5/input_chr3_subset.bam -f BAM -n H3K4me3
macs2 callpeak -t /usr/local/share/data/assignment_5/H3K27ac_chr3_subset.bam -c /usr/local/share/data/assignment_5/input_chr3_subset.bam -f BAM -n H3K27ac
```

### b. Understanding the peaks

Macs2 calls the peaks by analyzing the enrichments of sequences in the
genome compared to the input. The algorithm uses different measures to
ensure the enrichments are statistically significant and make biological
sense based on the sequence length, the expected peaks (narrow or broad)
and other parameters. We will focus on the H3K4me3 mark to visualize how
are the peaks are identified.

    ## 1. Copy the H3K4me3 .narrowPeak file to your local computer

    #scp bkirsh@137.82.55.186:/home/bkirsh/assignment5/H3K4me3_peaks.narrowPeak /mnt/c/Users/bkirs/Documents

    ## 2. Open IGV with and load the hg38 reference 

    ## 3. Load the following files:
    ### a. H3K4me3 bigwig file that you created in section 1 
    ### b. Input  bigwig file that you created in section 1 
    ## Note: remember to autoscale each file track!

    ## 4. Go to this position: * chr3:44,608,952-44,670,849 *

    #?# Compare the input and H3K4me3 signal trakcs (a and b), are there regions that you consider to have an enriched signal in the H3K4me3 track compared to the input? If yes, how many?- 0.5 pt

    #ANSWER: Yes, two.

    #?# Why do you consider those regions to have an enrichment? Would you need to confirm that the enrichment is statistically significant? If so, what do you think would be a way to test for the significance (you can mention the tools used on this assignment or explain conceptually what would you do to confirm an enrichment)? - 1 pt

    #ANSWER: In the H3K4me3 signal track, there are two distinctly sharp peaks with noticeably higher amplitude than the rest of the region. In comparison, the input track contains consistent and uniformly distributed background signal, without anbormally high amplitude regions or clear peak separation. Although a visual inspection of the H3K4me3 signal track enables peak visualization, further statistical testing is required to confirm their significance. I would use Bioconductor's ChipQC package or MACS and deepTools themselves to assess signal-to-noise ratio based on metrics such as strand cross-correlation (SE), irreproducibility discovery rate (IDR), and fingerprint plots. These quality control procedures help with distinguishing the ChIP signal from the background signal and can be calculated using the tools mentioned above.

    ## 5. Load into the same session the following files:
    ### c. H3K4me3 normalized using the input bigWig file that was pre-computed and you used on section 2
    ### d. H3K4me3 .narrowPeak file that you created in section 3a
    ## Note: remember to autoscale each file track!

    #scp bkirsh@137.82.55.186:/usr/local/share/data/assignment_5/H3K4me3_norm.bw /mnt/c/Users/bkirs/Documents

    ## 6. Go to the same position as in step 4

    #?# Does the region/s you thought had an enrichment show differences between the input and the H3K4me3 mark in the bamCompare (file c)? - 0.5 pt

    #ANSWER: Yes. Once again, the two peaks in the bamCompare file can be easily visualized because they lie on the positive axis as compared to the remaining background signal, which is located on the negative axis.

    #?# Are these visually enriched regions you selected, found in file d (narrowPeak)? What does that mean? - 1 pt

    #ANSWER: Yes. This means that the peaks are true peaks (i.e., statistically significant).

    #?# Add a screenshot of your IGV session right after this chunk, using Rmarkdown syntax - 1pt

    #ANSWER: Done, please find the screenshot below.

\*\*\* SCREENSHOT OF IGV\*\*\*

![Visualization of H3K4me3 signals using IGV
(chr3:44,608,952-44,670,849).](/Users/bkirs/Documents/School/BMEG400E/Assignments/igv_screenshot2.png)

## c. Peak enrichments

For this assignment, we are working with 3 different epigenetic marks:
H3K4me3, H3K27me3 and H3K27ac. Each of these marks a different
transcription state (i.e., activation or repression). Thus, to better
visualize the differences in the called peaks for each of these
epigenetic marks you will create a heatmap plot (like the one you
created on part 2) using their called peaks files (.narrowPeak or
.broadPeak).

``` bash

### Create 3 heatmaps following the specifications you used on part 2 (one for the peaks called for each epigenetic mark, but containing data from all three tracks). Use the peak files of each of the epigenetic marks as reference files. Use ONLY the non-input normalized files: H3K27ac_norm.bw H3K27me3_norm.bw H3K4me3_norm.bw
#?# Write the commands you used to compute the matrices and create the heatmaps below: - 3 pt

computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K27me3_norm.bw /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw -R H3K27me3_peaks.broadPeak -a=1000 -b=1000 -o matrix_norm_H3K27me3.gz

computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K27me3_norm.bw /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw -R H3K4me3_peaks.narrowPeak -a=1000 -b=1000 -o matrix_norm_H3K4me3.gz

computeMatrix scale-regions -S /usr/local/share/data/assignment_5/H3K27me3_norm.bw /usr/local/share/data/assignment_5/H3K4me3_norm.bw /usr/local/share/data/assignment_5/H3K27ac_norm.bw -R H3K27ac_peaks.narrowPeak -a=1000 -b=1000 -o matrix_norm_H3K27ac.gz


plotHeatmap -m matrix_norm_H3K27me3.gz --hclust 3 -o matrix_norm_clust_H3K27me3.png --colorMap Blues
plotHeatmap -m matrix_norm_H3K4me3.gz --hclust 3 -o matrix_norm_clust_H3K4me3.png --colorMap Blues
plotHeatmap -m matrix_norm_H3K27ac.gz --hclust 3 -o matrix_norm_clust_H3K27ac.png --colorMap Blues

#?# Add screenshots of the 3 heatmaps you got using the epigenetic marks' peak files as reference files. Add them after this code chunk in the following order: H3K4me3, H3K27ac, H3K27me3 - 1.5 pt

scp bkirsh@137.82.55.186:/home/bkirsh/assignment5/matrix_norm_clust_H3K27me3.png /mnt/c/Users/bkirs/Documents
scp bkirsh@137.82.55.186:/home/bkirsh/assignment5/matrix_norm_clust_H3K4me3.png /mnt/c/Users/bkirs/Documents
scp bkirsh@137.82.55.186:/home/bkirsh/assignment5/matrix_norm_clust_H3K27ac.png /mnt/c/Users/bkirs/Documents
```

***H3K4me3 screenshot***

![Heatmap of the distribution of input-normalized epigenetic marks
(H3K27me3, H3K4me3, and H3K27ac) relative to the peak file for
H3K4me3.](/Users/bkirs/Documents/matrix_norm_clust_H3K4me3.png)

***H3K27ac screenshot***

![Heatmap of the distribution of input-normalized epigenetic marks
(H3K27me3, H3K4me3, and H3K27ac) relative to the peak file for
H3K27ac.](/Users/bkirs/Documents/matrix_norm_clust_H3K27ac.png)

***H3K27me3 screenshot***

![Heatmap of the distribution of input-normalized epigenetic marks
(H3K27me3, H3K4me3, and H3K27ac) relative to the peak file for
H3K27me3.](/Users/bkirs/Documents/matrix_norm_clust_H3K27me3.png)

    #?# Do you see an overlap between the peaks of different epigenetic marks? Which epigenetic marks? - 1 pt

    #ANSWER: Yes, there are slight overlaps between the peaks of all three epigenetic marks. However, the ones between H3K4me3 and H3K27ac are particularly noticeable, especially between the TSS and TES of genes. 

    #?# Why do you think these epigenetic marks overlap? - 1 pt

    #ANSWER: Functional combinations of H3K4me3 and H3K27ac are associated with multimark motifs such as active promoters, while H3K27me3 indicates repressed or poised promoters. Thus, we can expect to see frequent co-localization of H3K4me3 and H3K27ac marks. Furthermore, broad histone marks such as H3K27me3 do not generally share motifs with narrow marks such as H3K27ac and H3K4me3. 

    #?# Explain the pattern you see for the peaks of H3K27me3, do they look similar or different to the other two marks? Why do you think this happens? - 1pt

    #ANSWER: The H3K27me3 signals are different from the other two because they are more uniformly spread out across the gene body. In contrast, the peaks for H3K4me3 and H3K27ac reach their maxima between the TSS and the TES of the gene body. 

    #?# Why do you think the borders of the elements have such clearly-defined borders with respect to the ChIP signal? -1 pt

    #ANSWER: I believe this is the case because different histone modifications have distinctive preferences in genomic localizations, which dictate the differences in their impact on gene expression. For example, H3K4me3 is associated with transcriptionally active gene promoter regions, thus having clearly-defined borders with respect to the ChIP signal between the TSS and TES of genes. Repressed genes have a much higher density of nucleosomes and can be marked by H3K27me3, which explains why this mark appears throughout the gene body.
