#  - no alignment
#  - no global recentering
#  - unwrap ONLY the dimer
#  - wrap solvent separately
#  - wrap dimer as intact fragments into the unit cell center
# Requires: pbctools (as per your help page)

############################################################
# Input files
set psf_file "../step3_input_nowat.psf"
set dcd_file "ck_S20.dcd"

# >>> IMPORTANT: define BOTH chains/segments of your dimer here <<<
set DIMER_STRING "not (water)"

# Output
set OUT_DCD "ck_pbc.dcd"
############################################################

# Load
mol new $psf_file type psf waitfor all
mol addfile $dcd_file type dcd waitfor all
set mol [molinfo top]
set nframes [molinfo $mol get numframes]
if {$nframes < 1} {
    puts "ERROR: no frames."
    return
}

# pbctools
if {[catch {package require pbctools} emsg]} {
    puts "ERROR: pbctools not available: $emsg"
    return
}

# Sanity: make sure dimer selection exists
set test_dimer [atomselect $mol $DIMER_STRING]
if {[$test_dimer num] == 0} {
    puts "ERROR: DIMER_STRING matched 0 atoms: '$DIMER_STRING'"
    return
}
$test_dimer delete

# --- Step 1: fix broken residues/fragments in the dimer on the first frame ---
# Join by residues first; if you still see breaks, try 'fragment' (see note below).
puts "Step 1: pbc join res on frame 0 for DIMER ..."
pbc join res -first 0 -last 0 -sel $DIMER_STRING -ref "name CA"

# --- Step 2: unwrap ONLY the dimer across the whole trajectory ---
# This removes big jumps over time without touching solvent, keeping the dimer continuous.
puts "Step 2: pbc unwrap DIMER across all frames ..."
pbc unwrap -sel $DIMER_STRING -first 0 -last [expr {$nframes-1}]

# --- Step 3: wrap solvent/ions (everything that's NOT the dimer) ---
# Keep residues intact; wrap to the unit cell center (default).
puts "Step 3: pbc wrap solvent/ions ..."
pbc wrap -sel "not ($DIMER_STRING)" -compound residue -center unitcell

# --- Step 4: wrap the dimer as intact fragments into the central image ---
# Using -compound fragment keeps each protein chain (fragment) whole;
# we do NOT center on COM; we just place them into the central unit cell.
puts "Step 4: pbc wrap DIMER as fragments into the unit cell ..."
pbc wrap -sel $DIMER_STRING -compound fragment -center unitcell

# --- Step 5: write out ---
puts "Step 5: writing -> $OUT_DCD ..."
animate write dcd $OUT_DCD waitfor all

puts "Done. Wrote: $OUT_DCD"
puts "If any residues in the dimer still show long bonds, rerun Step 1 with:"
puts "   pbc join fragment -first 0 -last 0 -sel $DIMER_STRING"
puts "And then repeat Steps 2–5."

exit
