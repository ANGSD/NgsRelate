(TeX-add-style-hook "ngsrelate_draft_v3"
 (lambda ()
    (LaTeX-add-bibitems
     "Patterson06"
     "Leutenegger"
     "Albrechtsen09"
     "Moltke11"
     "Albrechtsen10")
    (LaTeX-add-labels
     "eq1"
     "tab:prob"
     "tab:cond_prob"
     "tab:epsilon"
     "eq:xx")
    (TeX-run-style-hooks
     "verbatim"
     "amsmath"
     "geometry"
     "latex2e"
     "art11"
     "article"
     "a4paper"
     "11pt")))

