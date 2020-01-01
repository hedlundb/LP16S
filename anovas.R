library('Hmisc')
library('abind')

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

physeq = prune_samples(x=snn_physeq, samples=rownames(anovadat))

#physeq = prune_taxa(taxa_names(sn_physeq)[rowSums(otu_table(sn_physeq) > 1) > 1], sn_physeq)
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

#physeq = tax_glom(physeq, "Genus")

## Apply a p value correction to a series of anovas by turning them into a cube
## and then applying the fdr correction over the third dimension.
# pcorrect.3d = function(anovas)
# {
    # anovas.cube = do.call('abind', c(anovas, along=3))
    # anovas.pcorrect = plyr::alply(array(
                                        # apply(anovas.cube,
                                              # 3, p.adjust, method='fdr'),
                                            
                                        # dim=dim(anovas.cube),
                                        # dimnames = dimnames(anovas.cube))
                                  # , 3, .dims = TRUE)

    # anovas.pcorrect.df = do.call('rbind', lapply(names(anovas),
        # function(x)
        # {
            # z = anovas.pcorrect[x][[1]]
            # if(!is.null(z))
            # {
                # rownames(z) = paste0(x, '|', rownames(z))
            # }
            # return(z)
        # }))
    # return(anovas.pcorrect.df)
# }

# anovadat = dat[dat$Group == 'New',]

# anovas.location = tapply(1:nrow(anovadat), anovadat$Location,
    # function(rows)
    # {
        # regrdat = anovadat[rows,]
        
        # anov = do.call('cbind', setNames(lapply(colnames(diversity),
            # function(index)
            # {
                # model = aov(get(index) ~ Time, data=regrdat)
                # tukey = TukeyHSD(model)

                # ps = setNames(summary(model)[[1]][,5], rownames(summary(model)[[1]]))
                # ps = ps[-length(ps)]
                # out = c(ps, tukey[[1]][,4])
                # return(out)
            # }), colnames(diversity)))
        
       # return(anov)
    # })
    
# anovas.time = tapply(1:nrow(anovadat), anovadat$Time,
# function(rows)
# {
    # regrdat = anovadat[rows,]
    
    # anov = do.call('cbind', setNames(lapply(colnames(diversity),
        # function(index)
        # {
            # model = aov(get(index) ~ Location, data=regrdat)
            # tukey = TukeyHSD(model)

            # ps = setNames(summary(model)[[1]][,5], rownames(summary(model)[[1]]))
            # ps = ps[-length(ps)]
            # out = c(ps, tukey[[1]][,4])
            # return(out)
        # }), colnames(diversity)))
    
   # return(anov)
# })

# anovas.time.location = do.call('cbind', setNames(lapply(colnames(diversity),
        # function(index)
        # {
            # model = aov(get(index) ~ Time + Location, data=anovadat)
            # tukey = TukeyHSD(model)
            
            # ps = setNames(summary(model)[[1]][,5], rownames(summary(model)[[1]]))
            # ps = ps[-length(ps)]
            # out = c(ps, do.call('rbind',tukey)[,4])
            # return(out)
        # }), colnames(diversity)))

# anovas.location.pcorrect.df = pcorrect.3d(anovas.location)
# anovas.time.pcorrect.df = pcorrect.3d(anovas.time)

# write.csv(anovas.location.pcorrect.df, 'anovas.location.pcorrect.csv')
# write.csv(anovas.time.pcorrect.df, 'anovas.time.pcorrect.csv')

anovadat.physeq = prune_samples(x=snn_physeq, samples=rownames(anovadat))
community.mat = phyloseq:::veganifyOTU(anovadat.physeq)
        
dists = list(
jaccard = vegan::vegdist(community.mat, distance='jaccard'),
bray = vegan::vegdist(community.mat, distance='bray'),
unifrac = phyloseq:::fastUniFrac(anovadat.physeq, weighted=FALSE, normalized=TRUE, parallel = FALSE),
wunifrac = phyloseq:::fastUniFrac(anovadat.physeq, weighted=TRUE, normalized=TRUE, parallel = FALSE))

pairwise.adonis = function(physeq,data,var)
{
    combinations = with(data, combn(as.character(unique(get(var))), 2))
    lists = apply(combinations, 2,
        function(groups)
        {
            subdata = data[data[,var] %in% groups,]
            samples = rownames(subdata)
            subphyseq = prune_samples(x=physeq, samples=samples)
            submat = phyloseq:::veganifyOTU(subphyseq)
            subdists = list(
            jaccard = vegan::vegdist(submat, distance='jaccard'),
            bray = vegan::vegdist(submat, distance='bray'),
            unifrac = phyloseq:::fastUniFrac(subphyseq, weighted=FALSE, normalized=TRUE, parallel = FALSE),
            wunifrac = phyloseq:::fastUniFrac(subphyseq, weighted=TRUE, normalized=TRUE, parallel = FALSE))
            
            adon = do.call('rbind',lapply(subdists,
                function(d)
                {
                    na.omit(adonis(d ~ get(var), data=subdata)[[1]][,6])
                }))
            colnames(adon) = var
            
            simp = summary(simper(phyloseq:::veganifyOTU(subphyseq), subdata[,var]))[[1]]
            simp$contribution = diff(c(0,simp$cumsum))
            
            return(list(adon = adon, simp = simp))
        })
        
        adons = do.call('cbind',lapply(lists, function(x) return(x[[1]])))
        colnames(adons) = apply(combinations, 2, paste0, collapse='-')
        
        simps = do.call('cbind',lapply(lists, function(x) return(setNames(x[[2]]$contribution, rownames(x[[2]]))[rownames(lists[[1]][[2]])])))
        colnames(simps) = apply(combinations, 2, paste0, collapse='-')
        
        out = list(adons = adons, simps = simps)
        
        return(out)
}

adonis.loc = tapply(1:nrow(anovadat), anovadat$Location,
function(rows)
{
    regrdat = anovadat[rows,]
    samples = rownames(regrdat)
    regr.physeq = prune_samples(x=physeq, samples=samples)
    out = pairwise.adonis(regr.physeq, data = regrdat, var = 'Time')
    adons = out[[1]]
    colnames(adons) = paste0(regrdat$Location[1], ':', colnames(adons))
    simps = out[[2]]
    colnames(simps) = paste0(regrdat$Location[1], ':', colnames(simps))
    out = list(adons = adons, simps = simps)
    return(out)
})

adonis.loc.adonis = t(do.call('cbind',lapply(adonis.loc, function(x) return(x[[1]]))))
adonis.loc.simper = do.call('cbind',lapply(adonis.loc, function(x) return(x[[2]][rownames(adonis.loc[[1]][[2]]),])))

adonis.time = tapply(1:nrow(anovadat), anovadat$Time,
function(rows)
{
    regrdat = anovadat[rows,]
    samples = rownames(regrdat)
    regr.physeq = prune_samples(x=physeq, samples=samples)
    out = pairwise.adonis(regr.physeq, data = regrdat, var = 'Location')
    adons = out[[1]]
    colnames(adons) = paste0(regrdat$Time[1], ':', colnames(adons))
    simps = out[[2]]
    colnames(simps) = paste0(regrdat$Time[1], ':', colnames(simps))
    out = list(adons = adons, simps = simps)
    return(out)
})

adonis.time.adonis = t(do.call('cbind',lapply(adonis.time, function(x) return(x[[1]]))))
adonis.time.simper = do.call('cbind',lapply(adonis.time, function(x) return(x[[2]][rownames(adonis.time[[1]][[2]]),])))


adonis.loctime.adonis = rbind(adonis.time.adonis, adonis.loc.adonis)
adonis.loctime.simper = cbind(adonis.time.simper, adonis.loc.simper[rownames(adonis.time.simper),])

adonis.pw.df = do.call('rbind', lapply(dists,
    function(d)
    {
        z = adonis2(d ~ Time * Location + Individual, data = anovadat)
        return(setNames(unlist(z[,5])[1:4], c('Time', 'Location', 'Individual', 'Time&Location')))
    }))

write.csv(adonis.loctime.adonis, 'adonis.loctime.adonis.csv')
write.csv(adonis.loctime.simper, 'adonis.loctime.simper.csv')
write.csv(adonis.pw.df, 'adonis.pw.df.csv')