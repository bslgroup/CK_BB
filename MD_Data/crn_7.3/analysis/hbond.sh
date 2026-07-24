#!/bin/bash

replicates=(1 2 3)

# Format (UPDATED): AcceptorSegname AcceptorAtom DonorSegname DonorResid DonorAtom
pairs=(
  "HETA O2' PROA 194 ND1"
  "HETA O2' PROA 194 NE2"
  "HETB O2' PROB 194 ND1"
  "HETB O2' PROB 194 NE2"
)

# Build unique labels once (sanitize quotes in atom names for filenames)
for i in "${!pairs[@]}"; do
  pair=(${pairs[$i]})
  acceptor_seg=${pair[0]}
  acceptor_atom=${pair[1]}
  donor_seg=${pair[2]}
  donor_resid=${pair[3]}
  donor_atom=${pair[4]}

  # sanitize apostrophe in O2' -> O2p for filenames
  acceptor_atom_safe=${acceptor_atom//\'/p}

  pair_labels[$i]="dist_${acceptor_seg}${acceptor_atom_safe}-${donor_seg}${donor_resid}${donor_atom}"
done

for rep in "${replicates[@]}"; do
  psf_file="../$rep/step3_input_nowat.psf"
  pdb_file="../$rep/step3_input_nowat.pdb"
  dcd_file="../$rep/dcd/ck_wrapped.dcd"

  mkdir -p "../$rep/hbond"

  for i in "${!pairs[@]}"; do
    label="${pair_labels[$i]}"
    pair=(${pairs[$i]})

    acceptor_seg=${pair[0]}
    acceptor_atom=${pair[1]}
    donor_seg=${pair[2]}
    donor_resid=${pair[3]}
    donor_atom=${pair[4]}

    out_file="../$rep/hbond/${label}.dat"

    cat > hbond_dist.tcl << EOF
mol new $psf_file
mol addfile $pdb_file
mol addfile $dcd_file type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

# Move ALL atoms so ATP moves with the protein alignment
set all [atomselect top "all"]

# Reference selection (frame 0)
set ref [atomselect top "protein and alpha" frame 0]

# Mobile selection (update frame in loop)
set mob [atomselect top "protein and alpha"]

set donor    [atomselect top "protein and segname $donor_seg and resid $donor_resid and name $donor_atom"]
set acceptor [atomselect top "segname $acceptor_seg and resname ATP and name $acceptor_atom"]

# Optional sanity check (recommended)
if {[\$donor num] != 1 || [\$acceptor num] != 1} {
  puts "ERROR: selection not 1 atom. donor=[\$donor num], acceptor=[\$acceptor num]"
  quit
}

set outfile [open "$out_file" w]
set nframes [molinfo top get numframes]

for {set frame 0} {\$frame < \$nframes} {incr frame} {
  \$all frame \$frame
  \$mob frame \$frame
  \$donor frame \$frame
  \$acceptor frame \$frame

  # Align current protein(alpha) to frame 0, apply transform to ALL atoms
  set M [measure fit \$mob \$ref]
  \$all move \$M

  # Distance between donor atom and acceptor atom
  set dist [veclength [vecsub [measure center \$donor] [measure center \$acceptor]]]
  puts \$outfile "\$frame \$dist"
}
close \$outfile
quit
EOF

    echo "Running replicate $rep, pair $i (Acceptor: $acceptor_seg:$acceptor_atom, Donor: $donor_seg:$donor_resid:$donor_atom)..."
    vmd -dispdev text -e hbond_dist.tcl
    rm hbond_dist.tcl
  done
done

# --------- Gnuplot Plotting Section ---------
for i in "${!pair_labels[@]}"; do
  label="${pair_labels[$i]}"

  gnuplot <<- EOF
    set terminal postscript eps size 4,3.2 solid color enhanced lw 3.0 "Times-Bold" 30
    set termoption enhanced
    set encoding iso_8859_1
    set key vertical maxrows 3
    set bmargin 4
    set output "${label}.eps"
    set xlabel "Time (ns)" font ",30"
    set ylabel "Distance (\305)" offset 1,0 font ",30"
    set style line 1 lt 2 lw 2 lc rgb "gray" dt (5,5)
    set key top right
    set key width 1
    set key font "Arial,12"
    set xrange [0:1000]
    set yrange [0:20]
    set xtics 200
    set ytics 5

    plot "../1/hbond/${label}.dat" u (\$1*0.12*20):2 w l lw 0.5 lc rgb "blue" title "rep1", \
         "../2/hbond/${label}.dat" u (\$1*0.12*20):2 w l lw 0.5 lc rgb "black" title "rep2", \
         "../3/hbond/${label}.dat" u (\$1*0.12*20):2 w l lw 0.5 lc rgb "red" title "rep3", \
         3.5 w l ls 1 notitle
EOF

  epstopdf "${label}.eps"
  rm "${label}.eps"
  echo "Generated plot: ${label}.pdf"
done
