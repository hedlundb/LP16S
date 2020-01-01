#library('phyloseq')
#library('ape')
#library('Biostrings')

## importPhyseqFromPlaintext: Import and join physeq objects from filenames.
importPhyseqFromQiime2 = function(otutable=NULL,samdata=NULL,taxtable=NULL,tree=NULL, seqs=NULL)
{
    ## Parse files using specialized functions
    ## Prioritize for OTU table completeness.
    if (!is.null(otutable))
    {
        otumat = readTSVasOTUs(otutable)
        phy_out = phyloseq::otu_table(otumat, taxa_are_rows = TRUE)
    } else {
        stop('There is no otu table. Are you sure about this? Correct the code your damn self if you are.')
    }
    if (!is.null(taxtable))
    {
        taxmat = readTSVasTAX(taxtable)
        phy_out = phyloseq::merge_phyloseq(phy_out,phyloseq::tax_table(taxmat))
    }
    if (!is.null(samdata))
    {
        samdf = readTSVasSAMPLEDATA(samdata)
        phy_out = phyloseq::merge_phyloseq(phy_out,phyloseq::sample_data(samdf))
    }
    if (!is.null(tree))
    {
        phy_out = phyloseq::merge_phyloseq(phy_out,phyloseq::phy_tree(readNWKasTREE(tree)))
    }
    if (!is.null(seqs))
    {
        phy_out = phyloseq::merge_phyloseq(phy_out,phyloseq::refseq(readFASTAasSEQS(seqs)))
    }
    
    ## Create phyloseq object
    #taxa_names(physeq) = paste0('feature',c(1:ntaxa(physeq)))
    physeq = phyloseq::prune_taxa(phyloseq::taxa_sums(phy_out) > 0, phy_out)
    #attributes(physeq) = list(otu_table = attributes(physeq)[[1]],
    #                          tax_table = attributes(physeq)[[2]],
    #                          sam_data = attributes(physeq)[[3]],
    #                          phy_tree = attributes(physeq)[[4]],
    #                          ref_seq = attributes(physeq)[[5]],
    #                          class = attributes(physeq)[[6]],
    #                          feature_names = taxa_names(OTU))
    return(physeq)
}

## Read in a Qiime2 TSV feature table as a matrix
readTSVasOTUs = function(otutable)
{
    OTU = data.table::fread(otutable, sep='\t', header=TRUE, data.table = FALSE)
    ids = OTU[,1]
    otus = OTU[,-1]
    rownames(otus) = ids
    return(otus)
    #OTU = as.matrix(read.table(otutable,sep='\t',row.names=1,header=FALSE))
    #txt = scan(otutable,what='character',sep='\n',n=2,quiet=TRUE)[2]
    #splitxt = unlist(strsplit(txt,'\t'))
    #samnames = splitxt[2:length(splitxt)]
    #colnames(OTU) = samnames
    
    #return(OTU)
}

## Read in Qiime2 compatible sample metadata as a dataframe
readTSVasSAMPLEDATA = function(samdata)
{
    SAM = read.table(samdata,sep='\t',row.names=1,header=FALSE)
    txt = scan(samdata,what='character',sep='\n',n=2,quiet=TRUE)[1]
    splitxt = unlist(strsplit(txt,'\t'))
    samnames = splitxt[2:length(splitxt)]
    colnames(SAM) = samnames
    
    return(SAM)
}

## Read in a Qiime2 taxonomy table as a 7-column taxonomy matrix
readTSVasTAX = function(taxtable, ranks = c('Domain','Phylum','Class','Order','Family','Genus','Species'))
{
    TAX = data.table::fread(taxtable, header=TRUE, sep='\t')
    
    tt = do.call('rbind',lapply(strsplit(TAX$Taxon, ';'),
        function(x) setNames(c(x,rep('', max(0,length(ranks) - length(x))))[1:length(ranks)], ranks)
        ))
    colnames(tt) = ranks
    rownames(tt) = unlist(TAX[,1])
    # tt = read.table(taxtable,sep=c('\t'),row.names=1,header=FALSE)
    # CONF = tt[ncol(tt)]
    # tt = tt[-ncol(tt)]
    # taxa_list = strsplit(as.character(tt[,1]),';')
    # taxa_len = max(vapply(taxa_list, function(x) return(length(x)), FUN.VALUE=integer(1)))
    # taxa_fixed = lapply(taxa_list, function(x) return(c(x,unlist(replicate(taxa_len-length(x),"unclassified")))))
    # taxa_7level = lapply(taxa_fixed, function(x) return(x[1:nlevel]))
    
    # TAX = matrix(unlist(taxa_7level), nrow = length(taxa_list), ncol = 7, byrow=TRUE)
    # rownames(TAX) = rownames(tt)
    # colnames(TAX) = c('Domain','Phylum','Class','Order','Family','Genus','Species')
    
    return(tt)
}

readNWKasTREE = function(tree)
{
    TREE = ape::multi2di(ape::read.tree(tree))
    
    return(TREE)
}

readFASTAasSEQS = function(seqs)
{
    SEQS = Biostrings::readDNAStringSet(seqs)
    return(SEQS)
}