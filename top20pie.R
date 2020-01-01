library('ggplot2')
library('phyloseq')
source('tax_glom_fast.R')
#source('run_this.R')

taxpal = read.csv('Genus.taxa.pal.csv', header=TRUE, row.names=1, stringsAsFactors=FALSE)

physeq = snp_physeq

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

## Toy with the taxonomy so that the names are a little more consistent.
newtax = gsub('^D_.__', '', unlist(tax_table(physeq)[,6]))
subsetted = with(sample_data(physeq), {Time %in% c(0,3) & (Group %in% c('Old','New'))})
physeqp = subset_samples(physeq, subsetted)
# ## Tax glom (faster, cause it's slow otherwise.)
# taxglom = do.call('rbind', tapply(1:nrow(otu_table(physeq)), newtax,
    # function(x) colSums(otu_table(physeq)[x,])))[,subsetted]

taxglomt = tax_glom_fast(physeqp, 'Genus')
taxglom = otu_table(tax_glom_fast(physeqp, 'Genus'))
taxa_names(taxglom) = tax_table(taxglomt)[,6]

   
pieglom = do.call('cbind', tapply(1:ncol(taxglom), with(sample_data(physeqp), {paste(Group, Time, Location)}),
    function(x) if(length(x) > 1) rowSums(taxglom[,x]) else taxglom[,x]))

## Correct by dividing by sample size.
pieglomfixed = apply(pieglom, 2, function(x) x/sum(x))

## Sort by abundance.
pieglomfixed = pieglomfixed[rev(order(rowSums(pieglomfixed))),]

pieglomhead = pieglomfixed[1:19,]
pieglomhead = rbind(pieglomhead, 1 - colSums(pieglomhead))
rownames(pieglomhead)[nrow(pieglomhead)] = 'ZZ'
rownames(pieglomhead)[grepl('unclassified', rownames(pieglomhead))] = 'ZY'


taxa = sort(rownames(pieglomhead))
taxa_sub = gsub('ZY','Unclassified', taxa)
taxa_sub = gsub('ZZ','<5%', taxa_sub)

pieglomheadpal = setNames(taxpal[rownames(pieglomhead),], rownames(pieglomhead))
pieglomheadpal["ZZ"] = '#BDBDBD'
pieglomheadpal["ZY"] = '#CDCDCD'


library(ggplot2)
# Barplot
pieplots = apply(pieglomhead, 2,
    function(dat)
    {
        df = data.frame(group=rownames(pieglomhead), value = dat)
        bp = ggplot(df, aes(x="", y=value, fill = group)) +
            geom_bar(width = 1, stat = "identity", color='black') +
            coord_polar("y", start=0) +
            theme_bw() +
            theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            legend.text=element_text(family='Arial', colour = 'black'),
            legend.title=element_text(family='Arial', colour = 'black',face='bold'),
            legend.background = element_blank(),
            legend.key.size = unit(0.6,"line"),
            legend.position = 'none',
            axis.ticks = element_blank(),
            axis.title.y = element_blank(), #element_text(family='serif',size=10, colour = 'black', face='bold'),
            axis.title.x = element_blank(),
            plot.background = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            ) +
            scale_fill_manual(values = pieglomheadpal, labels = taxa_sub)
            return(bp)
    })

dir.create('piecharts')
lapply(1:length(pieplots),
    function(n)
    {
        fname = paste0('piecharts/', names(pieplots)[n], '.svg')
        svglite::svglite(fname)
            print(pieplots[[n]])
        dev.off()
    })
    
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 

plotlegend = g_legend(pieplots[[1]] + theme(legend.position = 'right'))
svglite::svglite('piecharts/legend.svg')
    grid::grid.draw(plotlegend)
dev.off()
