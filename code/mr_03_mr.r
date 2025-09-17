setwd('~/easymr')
source('~/easymr/code/function.r')

# arg
args = commandArgs(T)
start = as.numeric(args[1]); end = as.numeric(args[2])
range = c(start:end)

pair = read.csv('info/mr_pair_torun.csv')
pair = pair%>%arrange(to)
prev_from = prev_to = 'NULL'
# range = 700:nrow(pair)

for (i in range){
  if (i%%100==0){
    print(sprintf('%s/%s', i, nrow(pair)))
  }
  from = pair[i, 'from']
  to = pair[i, 'to']
  file_from = sprintf('~/db/gwas/%s/clump/input/%s', strsplit(from, '_')[[1]][1], from)
  file_to = sprintf('~/db/gwas/%s/clean/%s.txt.gz', strsplit(to, '_')[[1]][1], to)
  file_clump = sprintf('~/db/gwas/%s/clump/output/%s.clumped', strsplit(from, '_')[[1]][1], from)
  file_out = sprintf('data/harmo_dat/%s_TO_%s.rdata', from, to)
  if (file.exists(file_out)){next}
  if (from != prev_from) {
    df1 = read.delim(file_from)
    snp = read.table(file_clump, sep='', header=1)[,3]
  }
  if (to != prev_to) {df2 = read.delim(file_to)}
  mr_res = get_mr_res(df1, df2, snp)
  save(mr_res, file=file_out)
  rm(mr_res)
  prev_from = from
  prev_to = to
}
