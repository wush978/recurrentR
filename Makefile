all: README.html

%.md : %.Rmd
	Rscript -e "library(knitr);knit('$<', '$@')"

%.html : %.md recurrentR.bib
	pandoc --mathjax -s --bibliography=recurrentR.bib $< -o $@

clean:
	rm README.html

build:
	-rm -rf /tmp/recurrentR
	git clone . /tmp/recurrentR
	cd /tmp/recurrentR/recurrentR && Rscript -e "library(Rcpp);compileAttributes('.')" && Rscript -e "library(roxygen2);roxygenize('.', roclets=c('rd', 'namespace'))" && R CMD build .

