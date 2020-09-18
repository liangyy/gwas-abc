f = read.table('OMIM-LD-block.tsv.gz.tmp', sep = '\t', header = T, stringsAsFactors = F)
f = f[ f$is_omim, ]
chr = unlist(lapply(strsplit(f$lead_var, '_'), function(x) {x[1]}))
pos = unlist(lapply(strsplit(f$lead_var, '_'), function(x) {x[2]}))
out = data.frame(chromosome = chr, position = pos, lead_var = f$lead_var, trait = f$trait)
gz = gzfile('OMIM-LD-block.tsv.gz', 'w')
write.table(out, gz, sep='\t', row=F, col=T, quo=F)
close(gz)
