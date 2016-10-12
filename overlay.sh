#!/bin/bash

# ../../h750/bin/overlay -i signif_"$1".root -o signif_"$1".pdf -r '/(.*)/\1' \
#   --regex.xtitle='/(.*)/\1/(pT_.*)/\1 [GeV]/pT_([^ ]*) (.*)/p_{T}^{\1} \2/yAbs_yy/|y_{#gamma#gamma}|/N_j_(.*)/N_{jets}^{p_{T}#geq\1}/cosTS_yy/|cos #theta_{#gamma#gamma}*|/Dphi_j_j/|#Delta#phi_{jj}|' \
#   --regex.ytitle='|.*|s/#sqrt{s+b}' \
#   --yrange="0:$2" --widths=4

overlay2 signif.root -o signif.pdf \
  --yrange="0:5" --widths=4 -m '0.1:0.05:0.13:0.05' \
  --ticks-left --val-fmt="3.1f" --marker-size=2.2 --marker-color=1 \
  --xlabel-size=1.4 --xtitle-size=1.5 --xtitle-offset=1.0 \
  --ylabel-size=1.4 --ytitle-size=1.4 --ytitle-offset=1.0 \
  -r 'y/^.*/s\/#sqrt{s+b}' 'ng' \
  'nx/^(pT|m)_.*/& [GeV]' \
  'x/pT_([^ ]*)/p_{T}^{\1}' \
  'x/yAbs_yy/|y_{#gamma#gamma}|' \
  'x/cosTS_yy/|cos #theta_{#gamma#gamma}*|' \
  'x/N_j_(.*)/N_{ jets}^{ #geq\1 GeV}' \
  'x/Dphi_j_j/|#Delta#phi_{jj}|' \
  'x/m_jj/m_{jj}' \
  'x/Dy_j_j/|#Deltay_{jj}|'
