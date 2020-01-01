## Required packages: data.table, phyloseq

tax_glom_fast = function (physeq, taxrank = rank_names(physeq)[1])#, rename.otus = FALSE)
## Agglomerate taxa using data.table for the actual agglomeration (much faster.)
{
    if (length(taxrank) > 1)
    {
        showWarning('More than one taxrank supplied. Using the first.')
        taxrank = taxrank[1]
    }

    if(is.character(taxrank))
    {
        if (taxrank %in% phyloseq::rank_names(physeq))
        {
            ti = which(phyloseq::rank_names(physeq) == taxrank)
        } else {
            stop('Invalid taxrank supplied.')
        }
    } else if(is.integer(taxrank))
    {
        if (taxrank <= length(phyloseq::rank_names(physeq)) & taxrank > 0)
        {
            ti = taxrank
        } else {
            stop('Invalid taxrank supplied.')
        }
    }

    ## Construct data.table object
    om = data.table::data.table(attributes(physeq)[[1]])
    tm = data.table::data.table(attributes(physeq)[[2]])
    names = phyloseq::sample_names(physeq)
    .tr = data.table::data.table(tr=tm[,get(taxrank)])
    .tl = data.table::data.table(tl=phyloseq::taxa_names(physeq))
    colnames(.tr) = '.tr'
    colnames(.tl) = '.tl'
    proc = data.table::data.table(om, .tr, .tl)
    tm$.tl = .tl 
    #rename.otus = TRUE
    #if (!rename.otus)
    #{
    ## Agglomerate by taking the first OTUname of each taxon, then sum
    ## columns by each agglomerated otu name.
    proc[, .glom := .tl[1], by=.tr]
    summed = proc[, lapply((mget(names)),sum), by=".glom"]
    data.table::setkey(summed, '.glom')
    data.table::setkey(tm, '.tl')
    
    ## Construct output
    tout = as.matrix(tm[.tl %in% summed[,.glom],])
    rownames(tout) = tout[,'.tl']
    tout = tout[,-which(colnames(tout) == '.tl')]
    
    tout[,1:ncol(tout) > ti] = NA
    
    oout = as.matrix(summed[,-1])
    dimnames(oout) = list(unlist(summed[,1]),names)
    
    phyout = phyloseq::merge_phyloseq(phyloseq::otu_table(oout, taxa_are_rows=TRUE),
                                      phyloseq::tax_table(tout),
                                      phyloseq::refseq(physeq, errorIfNULL=FALSE),
                                      phyloseq::phy_tree(physeq, errorIfNULL=FALSE),
                                      phyloseq::sample_data(physeq, errorIfNULL=FALSE))
    
    return(phyout)
}

subset_then_tax_glom = function(physeq, rank)
{
    phyloseq::tax_glom(phyloseq::rarefy_even_depth(physeq, 100), rank)
}

subset_then_tax_glom_fast = function(physeq, rank)
{
    tax_glom_fast(phyloseq::rarefy_even_depth(physeq, 100), rank)
}

## Don't run these because tax_glom is very slow.
# microbenchmark::microbenchmark(tax_glom_fast(physeq, 'Phylum'), tax_glom_fast2(physeq, 'Phylum'), times=10L)
# microbenchmark::microbenchmark(subset_then_tax_glom (physeq, 'Genus'), subset_then_tax_glom_fast(physeq, 'Genus'), times=2L)

                                       # expr      min       lq     mean   median       uq      max neval cld
     # subset_then_tax_glom(physeq, "Phylum") 2459.992 2459.992 2600.651 2600.651 2741.309 2741.309     2   b
# subset_then_tax_glom_fast(physeq, "Phylum")  610.181  610.181  615.485  615.485  620.789  620.789     2   a
 
                                        # expr     min        lq       mean     median        uq       max neval cld
     # subset_then_tax_glom(physeq, "Genus") 25453.997 25453.997 25932.4705 25932.4705 26410.944 26410.944     2   b
# subset_then_tax_glom_fast(physeq, "Genus")   494.938   494.938   585.7255   585.7255   676.513   676.513     2   a
