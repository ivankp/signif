#!/bin/bash

../../h750/bin/overlay -i signif_"$1".root -o signif_"$1".pdf -r '/(.*)/\1' \
  --regex.xtitle='/(.*)/\1/(pT_.*)/\1 [GeV]/pT_([^ ]*) (.*)/p_{T}^{\1} \2/yAbs_yy/|y_{#gamma#gamma}|/N_j_(.*)/N_{jets}^{p_{T}#geq\1}/cosTS_yy/|cos #theta_{#gamma#gamma}*|/Dphi_j_j/|#Delta#phi_{jj}|' \
  --regex.ytitle='|.*|s/#sqrt{s+b}' \
  --yrange="0:$2" --widths=4
