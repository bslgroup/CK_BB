atomselect macro allprot {protein}
#all protein
atomselect macro allprota {protein and alpha}
#protein and alpha carbons
atomselect macro ecreg1 {segname PROA and resid 219 to 223 228 to 232 238 to 245 and name CA}
#all atoms in the intercellular region 1 with their alpha carbons
atomselect macro ecreg2 {segname PROB and resid 219 to 223 228 to 232 238 to 245 and name CA}
#all atoms in the intercellular region 2 with their alpha carbons


atomselect macro hel4 {protein and resid 213 to 259 and alpha}
atomselect macro hel5 {protein and resid 270 to 323 and alpha}
atomselect macro hel6 {protein and resid 328 to 369 and alpha}
atomselect macro hel10 {protein and resid 857 to 902 and alpha}
atomselect macro hel11 {protein and resid 913 to 965 and alpha}
atomselect macro hel12 {protein and resid 971 to 1013 and alpha}

