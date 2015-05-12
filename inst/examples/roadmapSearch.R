query <- roadmapSearch(sample = "H1 cell line")

query <- roadmapSearch(sample = "H1 cell line", experiment = "RRBS")

query <- roadmapSearch(accession = "GSM621357")

query <- roadmapSearch(sample = c("H1 cell line","H9 cell line"),
                       experiment = c("RRBS","H3K79me2"))
