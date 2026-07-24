gnuplot gnuplot-multi-rmsf.prm
#gnuplot gnuplot-mscl-angle-all2.prm

ps2epsi rmsf.ps
epstopdf rmsf.epsi
#ps2epsi mscl-angle-all2.ps
#epstopdf mscl-angle-all2.epsi

rm *ps *psi
#okular mscl-angle-all.pdf  mscl-angle-all2.pdf
