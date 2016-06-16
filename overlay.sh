#!/bin/bash

../../h750/bin/overlay -i sig.root -o sig.pdf -r '/(.*)/\1' \
  --regex.xtitle='/(.*)/\1/(pT_.*)/\1 [GeV]/pT_([^ ]*) (.*)/p_{T}^{\1} \2/yAbs_yy/|#Deltay_{#gamma#gamma}|/N_j_(.*)/N_{jets}^{p_{T}#geq\1}/cosTS_yy/|cos #theta_{#gamma#gamma}*|/Dphi_j_j/|#Delta#phi_{jj}|' \
  --regex.ytitle='|.*|s/#sqrt{s+b}' \
  --yrange='0:0.45'
