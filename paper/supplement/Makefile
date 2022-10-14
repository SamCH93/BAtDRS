all: pdf clean

# names of .Rnw .bib. and .bst files
FILE=batdrs-supplement

## generate .pdf file from .tex file
pdf: $(FILE).tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## clean
clean:  
	rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).bbl \
	$(FILE).out  $(FILE).brf  $(FILE).synctex.gz \
	rm -r ./figure/
