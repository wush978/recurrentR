#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
Rscript -e "library(knitr);knit('recurrentR-guide.Rnw', 'recurrentR-guide.tex')"
pdflatex recurrentR-guide.tex
./run_bibtex.sh
pdflatex recurrentR-guide.tex
pdflatex recurrentR-guide.tex
pdflatex recurrentR-guide.tex
rm *.aux *.bbl *.blg *.log *.out *.toc

# mv -f *.pdf ../inst/doc/
# cp -f *.Rnw ../inst/doc/
