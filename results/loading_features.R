print("Loading features...")
prefix = "riboprot"
# prefix = "translfact"
features_filename = paste0("features_", prefix, ".rds")
if (!file.exists(features_filename)) {
  # ribo_genes = read.csv2("~/projects/probrain/data/ribo_genes.csv", sep="\t")
  # ribo_genes[ribo_genes$grps!="translation_factors",]$entrez
  # ribo_genes[ribo_genes$grps=="translation_factors",]$entrez
  # head(ribo_genes)
  # dim(ribo_genes)
  # get entrez
  # entrez = read.csv2("riboprot_entrez.csv", sep="\t")[,1]
  entrez = read.csv2("translfact_entrez.csv", sep="\t")[,1]
  genes_to_update = entrez
  url = "http://epimed.univ-grenoble-alpes.fr/database/query/jobid"
  jobid = jsonlite::fromJSON(url)
  url = "http://epimed.univ-grenoble-alpes.fr/database/query/genes/update"
  body = list(jobid=jobid, symbols=paste(genes_to_update, collapse=", "), taxid=9606)
  response = httr::POST(url, body = body, encode = "form")

  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/query/jobstatus?jobid=", jobid)
  job = jsonlite::fromJSON(url)
  WAIT = TRUE
  while (job$status != "success" & WAIT) {
    print(paste0("job:", job$jobid, " ", job$status, " ", job$current, "/", job$total))
    system("sleep 2")
    job = jsonlite::fromJSON(url)
  }
  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/query/jobs?jobid=", jobid)
  foo = read.csv2(url, header=TRUE, sep=";", stringsAsFactors=FALSE)
  
  # get entrez
  entrez = foo$entrez
  url = "http://epimed.univ-grenoble-alpes.fr/database/query/jobid"
  jobid = jsonlite::fromJSON(url)
  url = "http://epimed.univ-grenoble-alpes.fr/database/query/genes/position"
  assembly = "GRCh38"
  format = "bed"
  taxid = 9606
  WAIT = TRUE
  body=list(jobid=jobid, symbols=paste(entrez, collapse=", "), assembly=assembly, positionType="unique",  taxid=taxid)
  response = httr::POST(url, body = body, encode = "form")
  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/query/jobstatus?jobid=", jobid)
  job = jsonlite::fromJSON(url)
  while (job$status != "success" & WAIT) {
    print(paste0("job:", job$jobid, " ", job$status, " ", job$current, "/", job$total))
    system("sleep 2")
    job = jsonlite::fromJSON(url)
  }
  url = paste0("http://epimed.univ-grenoble-alpes.fr/database/query/jobs?format=", format, "&jobid=", jobid)
  features = read.csv2(url, header=TRUE, sep=";", stringsAsFactors=FALSE)
  features[,4:5] = features[,5:4]
  colnames(features)[4:5] = colnames(features)[5:4]
  rownames(features) = features[,4]
  features[,1] = factor(features[,1], levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrMT"))
  features = features[order(features[,1], features[,2]),]
  saveRDS(features, features_filename)  
}
features = readRDS(features_filename)
