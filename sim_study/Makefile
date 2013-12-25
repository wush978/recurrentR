all: Wang2001.html Huang2004.html Huang2010.html

%.md : %.Rmd
	Rscript -e "require(knitr);knit('$<', '$@')"

%.html : %.md recurrentR.bib
	pandoc --mathjax -s --bibliography=recurrentR.bib $< -o $@

clean:
	rm README.html
