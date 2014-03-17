RNW = $(wildcard *.Rnw)
TEX = $(RNW:.Rnw=.tex)

all : recurrentR-include recurrentR-guide.pdf

recurrentR-guide.pdf : $(TEX)
	./build_pdf.sh

%.tex : %.Rnw
	Rscript -e "library(knitr);knit('$<')"

clean :
	-rm recurrentR-guide.pdf
	-rm $(TEX)
