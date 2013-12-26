all: README.html

%.md : %.Rmd
	Rscript -e "library(knitr);knit('$<', '$@')"

%.html : %.md recurrentR.bib
	pandoc --mathjax -s --bibliography=recurrentR.bib $< -o $@

clean:
	rm README.html
