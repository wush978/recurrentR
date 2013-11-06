all: README.html

README.md : README.Rmd
	Rscript -e "require(knitr);knit('$<', '$@')"

README.html : README.md recurrentR.bib
	pandoc --mathjax -s --bibliography=recurrentR.bib README.md -o $@

clean:
	rm README.html
