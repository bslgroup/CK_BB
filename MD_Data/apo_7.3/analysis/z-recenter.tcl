# recenter_to_origin.tcl
# Loads the PBC-fixed trajectory and shifts each frame so that the dimer COM is at {0 0 0}.
# Uses only core VMD cmds (measure center), no pbctools needed.

# ----------- INPUT / OUTPUT -------------
set psf_file "../step3_input_nowat.psf"
set dcd_in   "ck_pbc.dcd"     ;# output from the previous script
set dcd_out  "ck_wrapped.dcd"

# IMPORTANT: selection that includes BOTH chains of the dimer
# e.g.: "protein and (chain A B)" or "segname A0 B0"
# If the only protein is the dimer, "protein" is OK.
set DIMER_STRING "not (water or ion or segname HETE HETF)"
# ----------------------------------------

# Load
mol new $psf_file type psf waitfor all
mol addfile $dcd_in type dcd waitfor all
set mol [molinfo top]
set nframes [molinfo $mol get numframes]
if {$nframes < 1} {
    puts "ERROR: no frames found in $dcd_in"
    return
}

# Selections
set sel_all   [atomselect $mol "all"]
set sel_dimer [atomselect $mol $DIMER_STRING]
if {[$sel_dimer num] == 0} {
    puts "ERROR: DIMER_STRING matched 0 atoms: '$DIMER_STRING'"
    return
}

puts "Recentering dimer COM to origin for $nframes frames ..."
for {set f 0} {$f < $nframes} {incr f} {
    animate goto $f
    $sel_all frame $f
    $sel_dimer frame $f
    $sel_all update
    $sel_dimer update

    # COM of the dimer (compat-safe)
    set com [measure center $sel_dimer weight mass]
    # Shift whole system by -COM so dimer COM goes to {0 0 0}
    set shift [list \
        [expr {-[lindex $com 0]}] \
        [expr {-[lindex $com 1]}] \
        [expr {-[lindex $com 2]}] ]
    $sel_all moveby $shift
}

puts "Writing -> $dcd_out ..."
animate write dcd $dcd_out waitfor all
puts "Done. Output: $dcd_out"
puts "Note: This only translates coordinates; it doesn't re-wrap or alter the unit cell."

