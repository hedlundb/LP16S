library('phyloseq')
library('ggplot2')
library('picante')
library('Hmisc')
library('svglite')
library('abind')

physeq.faithpd = function(physeq, split = TRUE)
## Quick faith's PD calculation using picante and phyloseq.
## Pass it a phyloseq object with a tree.
{
    if (!split) {
        OTU = taxa_sums(physeq)
    }
    else if (split) {
        OTU = as(otu_table(physeq), "matrix")
        if (taxa_are_rows(physeq)) {
            OTU = t(OTU)
        }
    }
    
    picante::pd(samp = OTU, tree = phy_tree(physeq), include.root=FALSE)
}

# rarefaction.curves = function(physeq)
# {
    # OTU = as(otu_table(physeq), "matrix")
        # if (taxa_are_rows(physeq)) {
            # OTU = t(OTU)
        # }
        
    # rarecurve(OTU, step = 1, xlab = "Sample Size", ylab = "Species", label = TRUE)
# }

#Color vector for individual:
#pal.individual = c("8" = "#f6e8c3", "558" = "#bf812d", "881" = "#827159", "97" = "#4c2a02", "568" = "#8c6bb1", "633" = "#e7d4e8", "635" = "#762a83")
pal.individual = c("8" = "#f6e8c3", "558" = "#bf812d", "881" = "#827159", "97" = "#4c2a02", "568" = "#FF0000", "633" = "#009999", "635" = "#9FEE00")
#Color vector for time:
pal.time = c("0 mo." = "#c6dbef", "1 mo." = "#6baed6", "9 mos." = "#2171b5", "14 mos." = "#08306b", ">2 yr." = "#CCCCCC")
#Color vector for sampling location:
pal.location = c("oral cavity" = "#7fdbcc", "gills" = "#f9f695", "carapace" = "#67af6d", "cloaca" = "#cc9368", "FloridaSubstrate" = "#fa9fb5", "SharkReefSubstrate" = "#f768a1")

estimate_richness_faithpd = function(physeq, split = TRUE)
## Access faith's pd function to return a dataframe with all diversity metrics.
{
    richnesses = estimate_richness(physeq, split)
     
    faithpd = matrix(physeq.faithpd(physeq, split)[,1], ncol=1, dimnames=list(sample_names(physeq), 'Faiths.PD'))
    
    out = cbind(richnesses, faithpd)
    return(out)
}

diversity = estimate_richness_faithpd(snr_physeq)
diversity = read.csv('diversity.csv', row.names=1, header=TRUE)
ignore = sample_data(snr_physeq)$Location %in% c("FloridaSubstrate", "SharkReefSubstrate", "lesion")

dat = cbind(sample_data(snr_physeq), diversity)
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
#ignore = sample_data(snr_physeq)$Location %in% c("FloridaSubstrate", "SharkReefSubstrate", "lesion") | sample_data(snr_physeq)$Group == 'New'

#for(metric in colnames(diversity))
#{
    metrics = c('Observed', 'Simpson', 'Shannon', 'Faiths.PD')
    #message(metric)

    # gdat = dat

    # grid_for_oldpoints = expand.grid(x=unique(dat$Time), a=unique(dat$Location))

    # nrows = apply(grid_for_oldpoints, 1, function(x) return(nrow(dat[dat$Time == x[1] & dat$Location == x[2] & dat$Group == 'Old',])))

    # grid_for_removal = grid_for_oldpoints[nrows < 2,]

    #gdat$boxplot_metric = gdat[,metric]

    # for(i in 1:nrow(grid_for_removal))
    # {
        # if (!any((gdat$Time == grid_for_removal[i,1]) & (gdat$Location == grid_for_removal[i,2]) & (gdat$Group == 'Old')))
        # {
            # gdat = rbind(gdat, gdat[1,])
            # gdat[nrow(gdat), c('Individual', 'Observed', 'Chao1', 'se.chao1', 'ACE', 'se.ACE', 'Shannon', 'Simpson', 'InvSimpson', 'Fisher', 'Faiths.PD')] = NA
            # rownames(gdat)[nrow(dat)] = paste0('APPENDED_FAKEGROUP_',i)
            # gdat[nrow(gdat), "Time"] = grid_for_removal[i,1]
            # gdat[nrow(gdat), "Location"] = grid_for_removal[i,2]
            # gdat[nrow(gdat), "Group"] = "Old"
        # }
        # gdat[(gdat$Time == grid_for_removal[i,1]) & (gdat$Location == grid_for_removal[i,2]) & (gdat$Group == 'Old'), "boxplot_metric"] = NA
    # }

    #dat[dat$Group == 'Old',"Individual"] = dat[dat$Group == 'Old',"Individual"][1]
    # levels(dat$Individual)[levels(dat$Individual) == dat[dat$Group == 'Old',"Individual"][1]] = 'Old'

    confints = do.call('rbind',lapply(metrics,
        function(metric)
        {
            do.call('rbind', tapply(1:nrow(dat), dat$Time,
              function(x) {z = dat[x,]; do.call('rbind', tapply(1:nrow(z), z$Location,
                         function(xx)
                             {
                             zz = z[xx,];
                             zz$mean = mean(zz[,metric], na.rm=TRUE)
                             zz$sd = sd(zz[,metric], na.rm=TRUE)
                             zz$se = zz$sd/sqrt(nrow(zz))
                             if(length(na.omit(zz$se)) < 3)
                             {
                                zz$mean = NA
                             }
                             zz$ci.min = zz$mean - 2*zz$se 
                             zz$ci.max = zz$mean + 2*zz$se
                             zz$metric = metric
                             zz[1,]
                             })
                         )}))
        }))
    mgdata = do.call('rbind', lapply(metrics,
        function(metric)
        {
            z = dat
            z$boxplot_metric = z[,metric]
            z$metric = factor(metric, levels = metrics)
            return(z)
        }))
    g = ggplot(mgdata) +
    geom_boxplot(aes(x = factor(Time), y = boxplot_metric, fill = factor(Time)), width = 0.75) +
    #geom_errorbar(data = confints, aes(x = factor(Time), ymin=ci.min, ymax=ci.max), color="black", width=.25) +
    #geom_point(data = confints, aes(x = factor(Time), y=mean, fill = factor(Time)), color="black", shape=21, size = 0.5) +
    geom_point(aes(factor(Time), boxplot_metric, color = Individual), shape=20, size = 1, alpha=0.8) +
    geom_line(aes(factor(Time), boxplot_metric, group = Individual, color = LineInd),size=0.5, alpha=0.8) +
    
    facet_grid(factor(metric, levels = metrics)~factor(Location), scales='free', switch = 'y') +
    theme_bw() +
    coord_cartesian(clip='off') + 
    theme(text=element_text(family='Arial',size=8, color = 'black'),
          axis.text = element_text(family='Arial',size=8, color = 'black'),
          panel.spacing.x=unit(0.1, "lines"),
          panel.spacing.y=unit(0.3, "lines"),
          axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, family='Arial', colour = 'black'),
          legend.text=element_text(family='Arial',size=8, colour = 'black'),
          legend.title=element_text(family='Arial',size=8, colour = 'black',face='bold'),
          strip.text=element_text(family='Arial',size=8, colour = 'black', face='bold'),
          strip.background.x = element_blank(), #element_rect(colour='black', fill='white'),
          strip.background.y = element_blank(), #element_rect(colour='black', fill='white'),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(linetype=3, colour='lightgray'),
          panel.grid.minor.y = element_blank(),
          legend.position='none',
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(family='Arial',size=8, colour = 'black', face='bold')) +
    scale_color_manual(breaks=names(pal.individual),values=pal.individual) +
    scale_fill_manual(breaks=names(pal.time),values=pal.time) +
    xlab('Time') +
    ylab(NULL) +
    scale_y_continuous(limits = c(0, NA)) +
    guides(color = guide_legend(override.aes = list(linetype = 0)))
    
    #png(paste0(metric,'-boxplot.png'))
    svglite('combined-95ci.svg', width=(8.5)/2.54, height=(8.5)/2.54)
        print(g)
    dev.off()
#}

#write.csv(diversity,'diversity.csv')
