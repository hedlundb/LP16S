library('Hmisc')
library('DESeq2')

ignore = sample_data(sn_physeq)$Location %in% c("FloridaSubstrate", "SharkReefSubstrate", "lesion")

dat = cbind(sample_data(sn_physeq))
dat$Location = factor(capitalize(as.character(dat$Location)))
dat$Time = factor(c('0 mo.', '1 mo.', '9 mos.', '14 mos.', '>2 yr.')[4*(dat$Group == 'Old') + (dat$Time+1)], levels = c('0 mo.', '1 mo.', '9 mos.', '14 mos.', '>2 yr.'))

if (!is.null(ignore))
## Remove samples if we so choose. 
{
    dat = dat[!ignore,]
}

dat$LineInd = dat$Individual
dat$LineInd[dat$Group == 'Old'] = NA
dat$Time[dat$Group == 'Old'] = '>2 yr.'

anovadat = dat[dat$Group == 'New',]

physeq = prune_taxa(taxa_names(sn_physeq)[rowSums(otu_table(sn_physeq) > 1) > 1], sn_physeq)
ranks = colnames(tax_table(physeq))
## Make substitutions in the taxonomy.
tax_table(physeq)[,1] = gsub(paste0('^D_0__'), '',tax_table(physeq)[,1])
for(i in 2:length(ranks))
{
    tax_table(physeq)[,i] = gsub(paste0('^D_',i-1,'__'), '',tax_table(physeq)[,i])
    uncl = grepl('[Uu]nclassified',tax_table(physeq)[,i]) |
    grepl('[Uu]ncultured',tax_table(physeq)[,i]) |
    grepl('[Mm]etagenome',tax_table(physeq)[,i])
    
    tax_table(physeq)[uncl,i] = ''
    
    tax_table(physeq)[tax_table(physeq)[,i] == '',i] = gsub('Unclassified Unclassified ', 'Unclassified ', paste0('Unclassified ', tax_table(physeq)[tax_table(physeq)[,i] == '', i-1]))
}

physeq = prune_samples(x=physeq, samples=rownames(anovadat))
#physeq = tax_glom(physeq, "Genus")

fitdds = function(design)
{
    dedata = anovadat
    dedata$Time.Location = paste0(dedata$Time,'.', dedata$Location)
    dds = DESeqDataSetFromMatrix(countData = t(phyloseq:::veganifyOTU(physeq)), colData = dedata, design = design)
    geomean = function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    geomeans = apply(counts(dds), 1, geomean)
    dds = estimateSizeFactors(dds, geoMeans = geomeans)
    dds = DESeq(dds, test="Wald", fitType="parametric", modelMatrixType = 'expanded', betaPrior = TRUE)

    return(dds)
}

dds = fitdds(~ Individual + Location + Time)

## The variable "Time" is a time series; we look at each timepoint across the
## dataset (1 vs 2, 2 vs 3, 3 vs 4) and then we look at initial vs final (1 vs 4)
time.inorder = paste0('Time', gsub(' ', '.', levels(anovadat$Time)[sort(unique(as.integer(anovadat$Time)))]))
time.sequence = cbind(time.inorder[c(1, length(time.inorder))], rbind(time.inorder[2:length(time.inorder)-1], time.inorder[2:length(time.inorder)]))
time.comparison.within = apply(time.sequence, 2, function(x){z = rep(0, length(resultsNames(dds))); z[resultsNames(dds) == x[1]] = 1; z[resultsNames(dds) == x[2]] = -1; z})
rownames(time.comparison.within) = resultsNames(dds)
colnames(time.comparison.within) = apply(time.sequence, 2, paste, collapse='-')

## Secondary to time is location; we want to know: at each location, which SVs
## are changing along with time? Therefore, we control each test for location.
location.newlevels = paste0('Location', gsub(' ', '.', levels(anovadat$Location)[sort(unique(as.integer(anovadat$Location)))]))
time.location.comparison.within = do.call('cbind', lapply(location.newlevels, function(x) {z=time.comparison.within; z[x,] = 1; colnames(z) = paste0(x, ':', colnames(z)); z}))

## Reverse it; look at each "between location" comparison at each timepoint.
location.comparison.within = apply(combn(location.newlevels, 2), 2, function(x){z = rep(0, length(resultsNames(dds))); z[resultsNames(dds) == x[1]] = 1; z[resultsNames(dds) == x[2]] = -1; z})
rownames(location.comparison.within) = resultsNames(dds)
colnames(location.comparison.within) = apply(combn(location.newlevels, 2), 2, paste, collapse='-')

## Repeat.
location.time.comparison.within = do.call('cbind', lapply(time.inorder, function(x) {z=location.comparison.within; z[x,] = 1; colnames(z) = paste0(x, ':', colnames(z)); z}))

## Combine the contrasts into a table which we will iterate over.
contrasts.of.interest = cbind(time.location.comparison.within, location.time.comparison.within)

## Run the contrasts of interest.
results.of.contrasts = apply(contrasts.of.interest, 2, function(x) as.data.frame(results(dds, contrast = as.vector(x), independentFiltering=TRUE, alpha=0.05, pAdjustMethod="fdr", parallel=TRUE, cooksCutoff = FALSE)))
results.with.taxonomy = lapply(results.of.contrasts, function(x)
    {
        x$Genus = tax_table(physeq)[rownames(x),"Genus"]
        x$color = read.csv('Genus.taxa.pal.csv',row.names=1, stringsAsFactors = FALSE, header=TRUE)[x$Genus,]
        x$color[is.na(x$color)] = '#BDBDBD'
        x$color = paste0('"',x$color,'"')
        return(x)
    })

## Function to access sample names by parsing resultNames from a DESeq object
parse.resultNames = function(dds, resultName)
{
    mm = as.data.frame(attributes(dds)$modelMatrix)
    m1 = strsplit(resultName,':')[[1]]
    mustfactor = m1[which(m1 %in% colnames(mm))]
    orfactor = strsplit(m1[which(!m1 %in% colnames(mm))],'-')[[1]]
    
    mustcolumn = mm[,mustfactor] > 0
    orcolumn = rowSums(mm[,orfactor]) > 0
    
    keeprows = mustcolumn & orcolumn
    
    return(rownames(mm)[keeprows])
}

## Write results to file.
dir.create('DESeq2s')
invisible(lapply(names(results.with.taxonomy),
    function(fpref)
    {
        print(fpref)
        samnames = parse.resultNames(dds, fpref)
        tempdat = anovadat[samnames,]
        tempdat = tempdat[order(tempdat$Time, tempdat$Location),]
        tab = as.data.frame(na.omit(results.with.taxonomy[fpref][[1]]))
        tab = cbind(tab,otu_table(snn_physeq)[rownames(tab),rownames(tempdat)])
        tab = tab[rev(order(1-tab[,"padj"])),]
        
        ## Generate a header comment for the table.
        if (all(tempdat$Time ==  tempdat$Time[1]))
        {
            begins = tapply(1:nrow(tempdat), tempdat$Location, function(x) {z = tempdat[x,]; rownames(z)[1]})
            middles = tapply(1:nrow(tempdat), tempdat$Location, function(x) {z = tempdat[x,]; rownames(z)[ceiling(nrow(z)/2)]})
            ends = tapply(1:nrow(tempdat), tempdat$Location, function(x) {z = tempdat[x,]; rownames(z)[nrow(z)]})
            begins = begins[!is.na(begins)]
            middles = middles[!is.na(middles)]
            ends = ends[!is.na(ends)]
            
            isbegin = rownames(tempdat) %in% begins
            ismiddle = rownames(tempdat) %in% middles
            isend = rownames(tempdat) %in% ends
            
            header.text = rep('----------', nrow(tempdat))
            header.text[ismiddle] = as.character(tempdat$Location)[ismiddle]
            header.text[isend] = paste0(header.text[isend], '|')
            header.text[isbegin] = paste0('|', header.text[isbegin])
            
            header.comment = paste0("#\t\t\t\t\t\t\t\t\t", paste(header.text, collapse='\t'),sep='')
        }
        if (all(tempdat$Location ==  tempdat$Location[1]))
        {
            begins = tapply(1:nrow(tempdat), tempdat$Time, function(x) {z = tempdat[x,]; rownames(z)[1]})
            middles = tapply(1:nrow(tempdat), tempdat$Time, function(x) {z = tempdat[x,]; rownames(z)[ceiling(nrow(z)/2)]})
            ends = tapply(1:nrow(tempdat), tempdat$Time, function(x) {z = tempdat[x,]; rownames(z)[nrow(z)]})
            begins = begins[!is.na(begins)]
            middles = middles[!is.na(middles)]
            ends = ends[!is.na(ends)]
            
            isbegin = rownames(tempdat) %in% begins
            ismiddle = rownames(tempdat) %in% middles
            isend = rownames(tempdat) %in% ends
            
            header.text = rep('----------', nrow(tempdat))
            header.text[ismiddle] = as.character(tempdat$Time)[ismiddle]
            header.text[isend] = paste0(header.text[isend], '|')
            header.text[isbegin] = paste0('|', header.text[isbegin])
            
            header.comment = paste0("#\t\t\t\t\t\t\t\t\t", paste(header.text, collapse='\t'),sep='')
        }

        ## Write to file.
        fname1 = paste0('DESeq2s/', fpref, '.all.tsv')
        #cat(header.comment, '\n\t', file = fname1)
        #write.table(tab, fname1, sep='\t', quote=FALSE, append=TRUE)
                
        keeps = tab[,"padj"] < 0.05 & rowSums(tab[,samnames]) > 0
        if(any(keeps))
        {
            tab2 = tab[keeps,]
        } else {
            tab2 = matrix(nrow=0, ncol=ncol(tab),dimnames=list(NULL,colnames(tab)))
        }
        fname2 = paste0('DESeq2s/', fpref, '.filter.tsv')
        cat(header.comment, '\n\t', file = fname2)
        write.table(tab2, fname2, sep='\t', quote=FALSE, append=TRUE)
    }))