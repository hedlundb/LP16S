library('ggplot2')
library('phyloseq')
library('grid')
#library('pBrackets')
library('svglite')
library('Hmisc')

source('tax_glom_fast.R')
source('run_this.R')
physeq = snp_physeq
ranks = colnames(tax_table(physeq))

pal = read.csv('palette.txt', header=FALSE, stringsAsFactors=FALSE)[,1]

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

g_legend = function(a.gplot)
{
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}


## Draw barplots for each rank.
for (current_rank in ranks)
{
    message(paste0('Agglomerating at the ', tolower(current_rank), ' level.'))
    physeq_glom = tax_glom_fast(physeq, taxrank=current_rank)
    message('\tDone agglomerating. Drawing plot.')
    
    ## Combine taxa at < 5 percent
    lt5p_taxa = taxa_names(physeq_glom)[taxa_sums(physeq_glom) < 0.05]
    physeq_glom_lt5p = merge_taxa(physeq_glom, lt5p_taxa)

    ## Replace the phylum name with "<5%"
    tax_table(physeq_glom_lt5p)[is.na(tax_table(physeq_glom_lt5p)[,current_rank]),current_rank] = '<5%'
    
    ## Melt phyloseq object.
    mdf = psmelt(physeq_glom_lt5p)
    
    taxa = levels(mdf[,current_rank])
    taxa_sub = gsub('Unclassified','ZY', taxa)
    taxa_sub = gsub('<5%','ZZ', taxa_sub)
    mdf[,current_rank] = factor(mdf[,current_rank],levels = levels(mdf[,current_rank])[order(taxa_sub)])
    
    taxa_pal = setNames(pal[1:length(taxa)], taxa)
    
    taxa_pal[sort(taxa_sub) == 'ZZ'] = '#BDBDBD'
    taxa_pal[grepl('ZY', sort(taxa_sub))] = '#CDCDCD'
    palcols = setNames(taxa_pal, levels(mdf[,current_rank]))
    
    write.csv(matrix(palcols, dimnames=list(names(palcols), c('color'))), paste0(current_rank,'.taxa.pal.csv'), col.names=NULL, quote=FALSE)


    ## Sort dataframe based on multiple variables.
    mdf_sort = mdf[with(mdf,order(Group, Individual, Location, Time, Sample, Abundance, get(current_rank))),]
    
    ## Change substrate so we can plot it.
    levels(mdf_sort$Individual) = unique(c(levels(mdf_sort$Individual), 'Substrate'))
    mdf_sort[mdf_sort$Group == 'Substrate', 'Individual'] = 'Substrate'
    
    ## Sample order.
    mdf_sort$Sample.Index = factor(mdf_sort$Sample, levels=mdf_sort$Sample[!duplicated(mdf_sort$Sample)])
    Facet.Size = setNames(table(cumsum(!duplicated(mdf_sort$Individual))),unique(mdf_sort$Individual))/length(unique(mdf_sort[,current_rank]))
    mdf_sort$Facet.Index = factor(mdf_sort$Individual, levels = names(Facet.Size))
    
    ## Time labels
    mdf_sort = within(mdf_sort,
    {
    
        xlab = paste0(capitalize(substr(Location,1, 4)),' ',c('0 mo.', '1 mo. ', '9 mos. ', '14 mos.', '>2 yr.')[4 * (Group == 'Old') + (Time + 1)])
        xlab[Facet.Index == 'Substrate'] = gsub('Subfstrate', '', Location[Facet.Index == 'Substrate'])
    })
    
    g = ggplot(mdf_sort[order(mdf_sort[,current_rank]),], aes_string(x='Sample.Index', y='Abundance', fill=current_rank)) +
        geom_bar(stat="identity", position='stack', color='black') +
        facet_wrap( ~ factor(Facet.Index), nrow=1, strip.position='bottom', scales='free_x') + 
        theme_bw() +
        scale_y_continuous(limits=c(0,1.001),expand = c(0, 1e-3)) +
        scale_x_discrete(breaks=mdf_sort$Sample.Index, labels = mdf_sort$xlab) +
        scale_fill_manual(values=palcols) + 
        theme(text = element_text(family = 'Arial', colour = 'black', size=8),
              axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, family='Arial', colour = 'black'),
              axis.text.y = element_text(family='Arial', colour = 'black'),
              legend.text=element_text(family='Arial', colour = 'black',  size=8),
              legend.title=element_text(family='Arial', colour = 'black',face='bold'),
              legend.background = element_blank(),
              legend.position = 'none',
              legend.key.size = unit(0.6,"line"),
              strip.text=element_text(family='Arial', colour = 'black', face='bold'),
              strip.background = element_blank(),
              strip.placement = "outside",
              axis.ticks = element_blank(),
              axis.title.y = element_blank(), #element_text(family='serif',size=10, colour = 'black', face='bold'),
              axis.title.x = element_blank(),
             # plot.margin = unit(c(0.1,0.1,0.02,0.1), "null"),
              plot.background = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank()) #+
              #guides(fill=guide_legend(nrow=ceiling(length(unique(mdf_sort[,current_rank]))/ 58)))
        
    ## Fix facet sizing.
    gt = ggplotGrob(g)
    # Get the column index in the gt layout corresponding to the panels.
    panelI = gt$layout$l[grepl("panel", gt$layout$name)]
    # Replace the default panel widths with relative widths.
    gt$widths[panelI] <- unit(Facet.Size, "null")
    # Add extra width between panels
    gt$widths[panelI + 1] = unit(0.03, "cm")
    
    for(i in which(grepl("strip-b", gt$layout$name))){
        gt$grobs[[i]]$layout$clip <- "off"
    }

    svglite(paste0(current_rank, '-barplot_all.svg'), width=(18)/2.54, height=(8)/2.54, pointsize=8)
        grid.draw(gt)
    dev.off()
    
    g_leg = g_legend(g + theme(legend.position = 'bottom'))
    svglite(paste0(current_rank, '-barplot_legend.svg'), width=(40)/2.54, height = (18)/2.54)
        grid.draw(g_leg)
    dev.off()

}