gnuplot gnuplot-multi-rmsd.prm
#gnuplot gnuplot-rmsd.prm

ps2epsi rmsd.ps
epstopdf rmsd.epsi
#ps2epsi mscl-angle-all2.ps
#epstopdf mscl-angle-all2.epsi

rm *ps *psi 
#okular mscl-angle-all.pdf  mscl-angle-al.pdf
