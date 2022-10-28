all: pdf clean

# names of .Rnw .bib. and .bst files
FILE=batdrs
BIB=bibliography
BIBSTYLE=apalikedoiurl

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw $(BIB).bib $(BIBSTYLE).bst
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## generate .pdf file from .tex (without first running the tex recipe)
pdf2:
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## clean
clean:  
	rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).tex  $(FILE).bbl \
	$(FILE).out  $(FILE).brf  $(FILE).synctex.gz \
	rm -r #./figure/

cleanfig:
	rm ./figure/
