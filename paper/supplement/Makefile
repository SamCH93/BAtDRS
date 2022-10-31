all: pdf clean

# names of .Rnw .bib. and .bst files
FILE=batdrs-supplement

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## clean
clean:  
	rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).bbl \
	$(FILE).out  $(FILE).brf  $(FILE).synctex.gz  $(FILE).tex \
#	rm -r ./figure/
