#!/bin/bash

replicates=(1 2 3)

# Format: AcceptorResid AcceptorAtom1 AcceptorAtom2 AcceptorAtom3 DonorResid DonorAtom1 DonorAtom2 DonorAtom3
pairs=(
  "329 OD1 OD2 -   202 OG  -  -"
  "202 OG  -   -   69  ND1 NE2 -"
  "329 OD1 OD2 -   69  ND1 NE2 -"
#  "311 OG  -   -   242 ND1 NE2 -"
#  "311 OG  -   -   441 ND1 NE2 -"
#  "179 OE1 OE2 -   248 NH1 NH2 NE"
#  "179 OE1 OE2 -   441 ND1 NE2 -"
#  "179 OE1 OE2 -   151 NH1 NH2 NE"
#  "179 OE1 OE2 -   147 NH1 NH2 NE"
#  "321 OD1 OD2 -   248 NH1 NH2 NE"
)

# Build unique labels once
pair_labels=()
for i in "${!pairs[@]}"; do
  pair=(${pairs[$i]})
  acceptor_resid=${pair[0]}
  donor_resid=${pair[4]}
  pair_labels[$i]="salt_${acceptor_resid}-${donor_resid}"
done

for rep in "${replicates[@]}"; do
  psf_file="../$rep/step3_input_nowat.psf"
  pdb_file="../$rep/step3_input_nowat.pdb"
  dcd_file="../$rep/dcd/ck_wrapped.dcd"
  mkdir -p "../$rep/bridges"

  for i in "${!pairs[@]}"; do
    label="${pair_labels[$i]}"
    pair=(${pairs[$i]})

    acceptor_resid=${pair[0]}
    acceptor_atom1=${pair[1]}
    acceptor_atom2=${pair[2]}
    acceptor_atom3=${pair[3]}
    donor_resid=${pair[4]}
    donor_atom1=${pair[5]}
    donor_atom2=${pair[6]}
    donor_atom3=${pair[7]}

    out_file="../$rep/bridges/${label}.dat"

    cat > salt_mindist.tcl << EOF
mol new $psf_file
mol addfile $pdb_file
mol addfile $dcd_file type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

set all [atomselect top "protein"]
set ref [atomselect top "protein and alpha" frame 0]
set com [atomselect top "protein and alpha"]
EOF

    [[ "$donor_atom1" != "-" ]] && echo "set D1 [atomselect top \"protein and resid $donor_resid and name $donor_atom1\"]" >> salt_mindist.tcl
    [[ "$donor_atom2" != "-" ]] && echo "set D2 [atomselect top \"protein and resid $donor_resid and name $donor_atom2\"]" >> salt_mindist.tcl
    [[ "$donor_atom3" != "-" ]] && echo "set D3 [atomselect top \"protein and resid $donor_resid and name $donor_atom3\"]" >> salt_mindist.tcl
    [[ "$acceptor_atom1" != "-" ]] && echo "set A1 [atomselect top \"protein and resid $acceptor_resid and name $acceptor_atom1\"]" >> salt_mindist.tcl
    [[ "$acceptor_atom2" != "-" ]] && echo "set A2 [atomselect top \"protein and resid $acceptor_resid and name $acceptor_atom2\"]" >> salt_mindist.tcl
    [[ "$acceptor_atom3" != "-" ]] && echo "set A3 [atomselect top \"protein and resid $acceptor_resid and name $acceptor_atom3\"]" >> salt_mindist.tcl

    cat >> salt_mindist.tcl << EOF

set outfile [open "$out_file" w]
set num_steps [molinfo top get numframes]

for {set frame 0} {\$frame < \$num_steps} {incr frame} {
    \$all frame \$frame
EOF

    [[ "$donor_atom1" != "-" ]] && echo "    \$D1 frame \$frame" >> salt_mindist.tcl
    [[ "$donor_atom2" != "-" ]] && echo "    \$D2 frame \$frame" >> salt_mindist.tcl
    [[ "$donor_atom3" != "-" ]] && echo "    \$D3 frame \$frame" >> salt_mindist.tcl
    [[ "$acceptor_atom1" != "-" ]] && echo "    \$A1 frame \$frame" >> salt_mindist.tcl
    [[ "$acceptor_atom2" != "-" ]] && echo "    \$A2 frame \$frame" >> salt_mindist.tcl
    [[ "$acceptor_atom3" != "-" ]] && echo "    \$A3 frame \$frame" >> salt_mindist.tcl

    cat >> salt_mindist.tcl << EOF
    \$all move [measure fit \$com \$ref]
    set distances {}
EOF

    for d in D1 D2 D3; do
      for a in A1 A2 A3; do
        if [[ "${!d}" != "-" && "${!a}" != "-" ]]; then
          echo "    catch {lappend distances [veclength [vecsub [measure center \$$d] [measure center \$$a]]]}" >> salt_mindist.tcl
        fi
      done
    done

    cat >> salt_mindist.tcl << EOF
    if {[llength \$distances] > 0} {
        set mindist [lindex [lsort -real \$distances] 0]
        puts \$outfile "\$frame \$mindist"
    } else {
        puts \$outfile "\$frame NaN"
    }
}
close \$outfile
quit
EOF

    echo "Running replicate $rep, pair $i (Acceptor: $acceptor_resid, Donor: $donor_resid)..."
    vmd -dispdev text -e salt_mindist.tcl
    rm salt_mindist.tcl
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
    set key outside top right
    set key width 1
    set title "CRN 5.5"
    set key font "Arial,12"
    set xrange [0:1000]
    set yrange [0:20]
    set xtics 200
    set ytics 5

    plot "../1/bridges/${label}.dat" u (\$1*0.12*20):2 w l lw 0.5 lc rgb "blue" title "rep1", \
         "../2/bridges/${label}.dat" u (\$1*0.12*20):2 w l lw 0.5 lc rgb "black" title "rep2", \
         "../3/bridges/${label}.dat" u (\$1*0.12*20):2 w l lw 0.5 lc rgb "red" title "rep3", \
         4.5 w l ls 1 notitle
EOF

  epstopdf "${label}.eps"
  rm "${label}.eps"
  echo "Generated plot: ${label}.pdf"
done
