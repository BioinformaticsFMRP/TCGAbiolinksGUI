
query <- biOmics::biOmicsSearch("brain")

# Experiments
# "Microarray"-"ExpressionArray"-"ExonArray"-"RNASeq"
# "MiRNAMicroArray"-"Firehose"-"DNAMethylation "-"miRNASeq"-"RRBS"
# "ChipSeq"-"MRESeq"-"Rampage"-"DNAsequencing"
# "fiveC"-RepliSeq"-"Others"
query <- biOmics::biOmicsSearch("brain", experiment = "ExpressionArray")
