library('Hmisc')
library('stringi')
library('ggplot2')
library('phyloseq')
library('ape')
library('gridExtra')

source('/mnt/c/Ubuntu/Data/PROCEDURE/scripts/io.align.R')

## Function to eliminate vertical gaps in an alignment.
reduce.alignment = function(align)
{
    nts = do.call('rbind',strsplit(do.call('c', align), ''))
    allgap = unlist(apply(nts, 2, function(x) all(x == '-')))
    nts = nts[,!allgap]
    message(paste0('reduce.alignment removed ', sum(allgap), ' sites'))
    return(setNames(as.list(apply(nts, 1, paste0, collapse = '')), names(align)))
}

## Play with the input data so that we have fields with appealing names.
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

## Fix the taxonomy of the data so that it isn't a mess.
physeq = prune_taxa(taxa_names(snn_physeq)[rowSums(otu_table(sn_physeq) > 1) > 1], snn_physeq)
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
sample_data(physeq) = sample_data(anovadat)
#physeq = tax_glom(physeq, "Genus")


## Unpack DESeq2 files based on their names.
## Get significant SVs.
files = list.files('DESeq2s/')
deseqs = lapply(paste0('DESeq2s/', files), function(fname) data.frame(data.table::fread(cmd = paste0("grep -v ^\\# ", fname), header=TRUE, sep='\t', quote=''),row.names=1))
signifs = unique(unlist(lapply(deseqs, rownames)))
signifs = signifs[signifs %in% taxa_names(physeq)]

## Keep only significant SVs in the object.
signif_physeq = prune_taxa(signifs, physeq)

# fileinfo = do.call('rbind', strsplit(files, ':'))

# constantvar = vapply(grepl('^Time', fileinfo[,1]), function(x) if(x) 'Time' else 'Location', character(1))
# contrastvar = vapply(grepl('^Time', fileinfo[,1]), function(x) if(x) 'Location' else 'Time', character(1))

# contrast.untrimmed = do.call('rbind', strsplit(stringi::stri_replace_last_fixed(fileinfo[,2],'.filter.tsv', ''), '-'))

# contrast.trimmed = cbind(stringi::stri_replace_first_fixed(contrast.untrimmed[,1], contrastvar, ''), stringi::stri_replace_first_fixed(contrast.untrimmed[,2], contrastvar, ''))
# tapply(1:nrow(contrast.trimmed)
# setNames(levels(mdat[,"Time"]), gsub(' ', '.' ,levels(mdat[,"Time"])))

# colnames(fileinfo) = c('files', 'constantstring', 'comparestring', 'constant')

## Write significant rep seqs to file.
# seqs = setNames(as.list(t(as.data.frame(refseq(signif_physeq)))), names(refseq(signif_physeq)))
# write.align(seqs, 'DESeq16S.fna')

## Do mothur stuff here.
## (Instead use arb-silva.de? It's the same thing.)
#system('sh makeTreeFromDESeqs.sh')

## Read tree, make it pretty, and extract top-to-bottom tiporder.
tree = phytools::midpoint.root(ape::ladderize(ape::read.tree('DESeq16S.silva.filter.fasta.treefile')))
tiporder = rev(tree$edge[tree$edge[,2] <= Ntip(tree),2])
tiplabels = tree$tip.label[tiporder]

## Get phyloseq-melt object
mdat = phyloseq::psmelt(signif_physeq)
mdat$otutip = factor(mdat$OTU, levels = tiplabels)
mdt = data.table::as.data.table(mdat)
#avemdt = mdt[,mean(Abundance), by=.(OTU,otutip,Location,Time)]
avemdat = do.call('rbind', tapply(1:nrow(mdat), paste0(mdat$Location,mdat$Time),
                                  function(x) {z = mdat[x,]; do.call('rbind', tapply(1:nrow(z), z$OTU,
                                                                                     function(xx) {zz = z[xx,]; zz$mean = mean(zz$Abundance); zz[1,]}))}))
avemdat$log10cpm = log10(avemdat$mean)
avemdat$log10cpm[avemdat$log10cpm == -Inf] = 0

## more playing with the taxonomy to make it fit the page.
genlabels = stringi::stri_replace_first_fixed(tax_table(signif_physeq)[levels(avemdat$otutip),"Genus"], 'Unclassified', 'Unc.')
genlabels = stringi::stri_replace_first_fixed(genlabels, 'Candidatus', 'Ca.')
genlabels = stringi::stri_replace_first_regex(genlabels, '^bacterium ', '')
genlabels = rev(setNames(vapply(genlabels, function(x) {substr(x, start=1,stop=1) = toupper(substr(x, start=1,stop=1)); return(x)}, character(1)), rownames(genlabels)))

## Grab the family level
taxfam = rev(stringi::stri_replace_first_fixed(tax_table(signif_physeq)[levels(avemdat$otutip),"Family"], 'Unclassified', 'Unc.'))
family_colors = read.csv('Family.taxa.pal.csv',sep=',',header=TRUE, row.names=1,comment.char='')

famlabels = setNames(taxfam, names(genlabels))
famlabelindex = setNames(cumsum(c(TRUE,unlist(lapply(2:length(taxfam), function(x) taxfam[x] != taxfam[x-1])))), names(famlabels))

avemdat$famlabelindex = famlabelindex[avemdat$otutip]
avemdat$Location = factor(avemdat$Location)
g = ggplot(data = avemdat, aes(x = Time, y = otutip, fill = log10cpm), clip = 'off') +
    geom_tile(width=0.9, height=0.9) +
    scale_fill_gradientn(colors = c("#F3F3F3", "#FFB913", "#8B2323")) +
    scale_y_discrete(breaks = levels(avemdat$otutip), labels = genlabels) +
    facet_grid(famlabelindex ~ Location, scales = "free_y", space = "free_y", labeller = labeller(famlabelindex = setNames(famlabels,famlabelindex), Location = avemdat$Location)) +
    coord_cartesian(ylim = c(0,length(genlabels))) +
    theme_bw() +
    coord_cartesian(expand=FALSE, clip = "off") +
    theme(text=element_text(family='Arial',size=8, colour = 'black'),
        panel.spacing.x=unit(0.3, "lines"),
        axis.text = element_text(family='Arial',size=8, colour = 'black'),
        axis.text.x = element_text(angle = -90, size=8, hjust = 0, vjust=0.5, family='Arial', colour = 'black'),
        legend.text=element_text(family='Arial',size=8, colour = 'black'),
        legend.title=element_text(family='Arial',size=8, colour = 'black',face='bold'),
        strip.text=element_text(family='Arial',size=8, colour = 'black', face='bold'),
        strip.text.y=element_text(family='Arial',size=11, colour = 'black',angle=0,hjust=0,vjust=1),
        strip.background.x = element_blank(),#element_rect(colour='black', fill='white'),
        strip.background.y = element_rect(colour='transparent', fill='gray'),
        strip.placement = "outside",
        panel.spacing.y = unit(0.1, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),#element_line(linetype=3, colour='lightgray'),
        panel.grid.minor.y = element_blank(),
        legend.position='none',
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
        
        
gt = ggplotGrob(g)
for(i in which(grepl("strip", gt$layout$name)))
{
        gt$grobs[[i]]$layout$clip <- "off"
}

labelz = which(grepl("axis-l", gt$layout$name))
firstlabel = min(labelz)

strips = which(grepl("strip-r", gt$layout$name))
firststrip = min(strips)

# gt$layout$t[strips] = gt$layout$t[labelz]
# gt$layout$b[strips] = gt$layout$b[labelz]
gt$layout$z[strips] = (1:length(gt$layout$z[strips]))+1
for(i in 1:length(strips))
{
    grobsize = sum(famlabelindex == i)
    gt$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]$vjust = unit((1-(0.5/grobsize))-(0.1/(grobsize/2)), 'npc')
    gt$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]$gp$fontsize = 6
    mycolor = family_colors[as.character(famlabels[!duplicated(famlabelindex)][i]),"color"]
    if (!is.na(mycolor))
    {
        gt$grobs[[strips[i]]]$grobs[[1]]$children[[1]]$gp$fill = as.character(mycolor)
        gt$grobs[[strips[i]]]$grobs[[1]]$children[[1]]$gp$color = as.character(mycolor)
    } else {
         gt$grobs[[strips[i]]]$grobs[[1]]$children[[1]]$gp$fill = "#BDBDBD"
    }
}

svglite::svglite('time-per-location.heatmap.svg', width=(16)/2.54, height=12)
    grid.draw(gt)
dev.off()
dev.off()


## location per time

g = ggplot(data = avemdat, aes(x = Location, y = otutip, fill = log10cpm), clip = 'off') +
    geom_tile(width=0.9, height=0.9) +
    scale_fill_gradientn(colors = c("#F3F3F3", "#FFB913", "#8B2323")) +
    scale_y_discrete(breaks = levels(avemdat$otutip), labels = genlabels) +
    facet_grid(famlabelindex ~ Time, scales = "free_y", space = "free_y", labeller = labeller(famlabelindex = setNames(famlabels,famlabelindex), Time = avemdat$Time)) +
    coord_cartesian(ylim = c(0,length(genlabels))) +
    theme_bw() +
    coord_cartesian(expand=FALSE, clip = "off") +
    theme(text=element_text(family='Arial',size=8, colour = 'black'),
        panel.spacing.x=unit(0.3, "lines"),
        axis.text = element_text(family='Arial',size=8, colour = 'black'),
        axis.text.x = element_text(angle = -90, size=8, hjust = 0, vjust=0.5, family='Arial', colour = 'black'),
        legend.text=element_text(family='Arial',size=8, colour = 'black'),
        legend.title=element_text(family='Arial',size=8, colour = 'black',face='bold'),
        strip.text=element_text(family='Arial',size=8, colour = 'black', face='bold'),
        strip.text.y=element_text(family='Arial',size=11, colour = 'black',angle=0,hjust=0,vjust=1),
        strip.background.x = element_blank(),#element_rect(colour='black', fill='white'),
        strip.background.y = element_rect(colour='transparent', fill='gray'),
        strip.placement = "outside",
        panel.spacing.y = unit(0.1, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),#element_line(linetype=3, colour='lightgray'),
        panel.grid.minor.y = element_blank(),
        legend.position='none',
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
        
        
gt = ggplotGrob(g)
for(i in which(grepl("strip", gt$layout$name)))
{
        gt$grobs[[i]]$layout$clip <- "off"
}

labelz = which(grepl("axis-l", gt$layout$name))
firstlabel = min(labelz)

strips = which(grepl("strip-r", gt$layout$name))
firststrip = min(strips)

# gt$layout$t[strips] = gt$layout$t[labelz]
# gt$layout$b[strips] = gt$layout$b[labelz]
gt$layout$z[strips] = (1:length(gt$layout$z[strips]))+1
for(i in 1:length(strips))
{
    grobsize = sum(famlabelindex == i)
    gt$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]$vjust = unit((1-(0.5/grobsize))-(0.1/(grobsize/2)), 'npc')
    gt$grobs[[strips[i]]]$grobs[[1]]$children[[2]]$children[[1]]$gp$fontsize = 6
    mycolor = family_colors[as.character(famlabels[!duplicated(famlabelindex)][i]),"color"]
    if (!is.na(mycolor))
    {
        gt$grobs[[strips[i]]]$grobs[[1]]$children[[1]]$gp$fill = as.character(mycolor)
        gt$grobs[[strips[i]]]$grobs[[1]]$children[[1]]$gp$color = as.character(mycolor)
    } else {
         gt$grobs[[strips[i]]]$grobs[[1]]$children[[1]]$gp$fill = "#BDBDBD"
    }
}

svglite::svglite('location-per-time.heatmap.svg', width=(16)/2.54, height=12)
    grid.draw(gt)
dev.off()
dev.off()