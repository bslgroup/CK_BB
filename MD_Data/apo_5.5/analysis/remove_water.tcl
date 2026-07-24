# ============================================================
# remove_tip3_min.tcl
# Minimal psfgen script to remove all TIP3 waters (OH2 atom)
# and write new PSF/PDB files.
#
# Usage:
#   vmd -dispdev text -e remove_tip3_min.tcl
# ============================================================

# --- Define input and output files ---
set inpsf  "step3_input.psf"
set inpdb  "step3_input.pdb"
set outpsf "step3_input_nowat.psf"
set outpdb "step3_input_nowat.pdb"

# --- Load structure into psfgen and VMD ---
package require psfgen
resetpsf
readpsf  $inpsf
coordpdb $inpdb

# Load into VMD for atom selection
mol load psf $inpsf pdb $inpdb

# --- Select all TIP3 water oxygens ---
set badwater [atomselect top "resname TIP3 and name OH2"]

# --- Delete those residues from psfgen (one per water) ---
foreach segid [$badwater get segid] resid [$badwater get resid] {
    delatom $segid $resid
}

# --- Write new files without TIP3 waters ---
writepsf $outpsf
writepdb $outpdb

puts "\n✅ Wrote new files without TIP3 waters:"
puts "  $outpsf"
puts "  $outpdb"
quit
