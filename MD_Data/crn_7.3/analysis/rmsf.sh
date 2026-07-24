#!/bin/bash

for j in `seq 0 0`
# loops through 2 system directories (change this number if you have more or less directories for example you have 5 directories the value goes for 0 to 4)
do

cat > RMSF_Apro.tcl << EOF
#Writes a RMSF scrpits for a single chain/segment you can combine this scrpit if you have one segment/chain in the system

mol new step3_input_nowat.psf
mol addfile step3_input_nowat.pdb
mol addfile dcd/ck_wrapped.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

source macros_ck.tcl

set all [atomselect top "allprot"]
set sel0 [\$all num]
set sel [atomselect top "proa1"]
set ref [atomselect top "proa1" frame 0]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {\$frame < \$num_steps} {incr frame} {
    \$all frame \$frame
    \$sel frame \$frame

set trans [measure fit \$sel \$ref]
\$all move \$trans
}

for {set i 0} {\$i < [\$sel num]} {incr i} {
#set rmsf [measure rmsf \$sel first 0 last -1 step 5]
     set rmsf [measure rmsf \$sel first 0 last -1]
     puts stderr "[expr {\$i+1}] \t [lindex \$rmsf \$i]"
#prints overall RMSF for each residue in the chain/segment
}

quit
EOF

cat > RMSF_Bpro.tcl << EOF

mol new step3_input_nowat.psf
mol addfile step3_input_nowat.pdb
mol addfile dcd/ck_wrapped.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

source macros_ck.tcl

set all [atomselect top "allprot"]
set sel0 [\$all num]
set sel [atomselect top "prob1"]
set ref [atomselect top "prob1" frame 0]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {\$frame < \$num_steps} {incr frame} {
    \$all frame \$frame
    \$sel frame \$frame

set trans [measure fit \$sel \$ref]
\$all move \$trans
}

for {set i 0} {\$i < [\$sel num]} {incr i} {
#     set rmsf [measure rmsf \$sel first 0 last -1 step 5]
     set rmsf [measure rmsf \$sel first 0 last -1]
     puts stderr "[expr {\$i+1}] \t [lindex \$rmsf \$i]"
}

quit
EOF

vmd -dispdev text -e RMSF_Apro.tcl 2> rmsf_a.txt

vmd -dispdev text -e RMSF_Bpro.tcl 2> rmsf_b.txt

done
