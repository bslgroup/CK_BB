
mol new step3_input_nowat.psf
mol addfile step3_input_nowat.pdb
mol addfile dcd/ck_wrapped.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

source macros_ck.tcl

set all [atomselect top "allprot"]
set sel0 [$all num]
set sel [atomselect top "prob1"]
set ref [atomselect top "prob1" frame 0]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {$frame < $num_steps} {incr frame} {
    $all frame $frame
    $sel frame $frame

set trans [measure fit $sel $ref]
$all move $trans
}

for {set i 0} {$i < [$sel num]} {incr i} {
#     set rmsf [measure rmsf $sel first 0 last -1 step 5]
     set rmsf [measure rmsf $sel first 0 last -1]
     puts stderr "[expr {$i+1}] \t [lindex $rmsf $i]"
}

quit
