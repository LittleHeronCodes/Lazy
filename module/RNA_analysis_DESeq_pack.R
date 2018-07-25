### DESeq RNA seq analysis pack ###


codingGenes = grch38[which(grch38$biotype == 'protein_coding'),]

# only protein coding genes
pcgidx = which(gsub('.[0-9]+$', '', rownames(rawCountMat)) %in% codingGenes$ensgene)

pcgCountMat = rawCountMat[pcgidx,]
ttt = apply(pcgCountMat, 1, function(v) sum(v==0) >=5)
count.filt = pcgCountMat[!ttt,]


cds = DESeqDataSetFromMatrix(countData = count.filt, colData = design, design = ~ modelType)
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)
countNorm = counts(cds, normalized =T)	# CPM (counts per million... i think)




