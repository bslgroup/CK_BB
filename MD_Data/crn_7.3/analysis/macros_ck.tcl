
atomselect macro allprot {protein}
#all protein
atomselect macro allprota {protein and alpha and not resid 81 to 104}
 #all protein and alpha carbons

atomselect macro proa1 {protein and segname PROA and alpha}
 #protein indiviual chain
atomselect macro prob1 {protein and segname PROB and alpha}


atomselect macro proa2 {protein and (segname PROA and helix) and alpha}
 #helix of chain indiviual
atomselect macro prob2 {protein and (segname PROB and helix) and alpha}

