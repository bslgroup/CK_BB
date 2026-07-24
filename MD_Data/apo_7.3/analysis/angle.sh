#!/bin/bash
#Ahmed Shubbar - UARK
#this scrpits includes anlysis for inter protomer angles , requires macros tcl file

cat >EC-angle.tcl << EOF

mol new step3_input_nowat.psf
mol addfile step3_input_nowat.pdb
mol addfile dcd/ck_wrapped.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

source macros_angle.tcl

set all [atomselect top "allprot"]
set EC_REG1 [atomselect top "ecreg1"]
set EC_REG2 [atomselect top "ecreg2"]
set com [atomselect top "allprota"]
set ref [atomselect top "allprota" frame 0 ]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {\$frame < \$num_steps} {incr frame} {
    \$all frame \$frame
    \$EC_REG1 frame \$frame
    \$EC_REG2 frame \$frame

 # Align to reference
    \$all move [measure fit \$com \$ref]

    set h1h2 [expr 180*(acos([vecdot [lindex [measure inertia \$EC_REG1] 1 2] [lindex [measure inertia \$EC_REG2] 1 2]]))/acos(-1)]
#Calculate the angle between  EC_REG1 and EC_REG2, use the dot product of their inertia tensors (CoM) and converting the angle from radians to degrees
    if {\$h1h2 > 90} {set h1h2 [expr 180-\$h1h2]}
#This command checks if the angle between the two regions is greater than 90 degrees. 
#If it is, then the angle is subtracted from 180 degrees.
# This is because the angle between two vectors can be either positive or negative, 
    puts stderr "\$frame \$h1h2"
#This command prints the frame number and the angle between the two regions to the standard error stream
}

quit
EOF


vmd -dispdev text -e EC-angle.tcl 2> EC-angle.txt
#runs the EC angle scrpit from above generated
