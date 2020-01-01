#library('phyloseq')
#library('ggplot2')
#library('picante')
#library('vegan')

source('importTSVPhyseq.R')

## Function used to generate the seed at the time of execution.
## Don't uncomment this unless you want to re-do the random seed.
#set.seed(as.integer(runif(1)*2e9))
#rseed = .Random.seed
#seed_date = date()
#save(".Random.seed",file="random_seed.RData")

## For reproducibility
#load("random_seed.RData")
#rseed = .Random.seed
#set.seed(rseed)


## NOTE: Tree will need to be re-done after samples are filtered.
physeq = importPhyseqFromQiime2(otutable = 'exported/plaintext_feature_table.tsv',
                                samdata = 'sample-metadata.tsv',
                                taxtable = 'exported/silva-taxonomy.tsv',
                                tree = 'exported/tree.nwk',
                                seqs = 'exported/dna-sequences.fasta')
source('applyFilter.R')


sn_physeq = phyloseq::prune_taxa(phyloseq::taxa_sums(s_physeq) > 0, s_physeq) # Non-normalized or rarefied
snn_physeq = phyloseq::transform_sample_counts(sn_physeq, function(x) 1E6 * x/sum(x)) # Normalized
snp_physeq = phyloseq::transform_sample_counts(sn_physeq, function(x) x/sum(x)) # Percent
snr_physeq = phyloseq::rarefy_even_depth(sn_physeq) # Rarefied
