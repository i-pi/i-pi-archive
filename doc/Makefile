.PHONY : all clean
all : manual.pdf

%.pdf : %.tex %.bbl
	while (pdflatex $*.tex ; \
	grep -q "Rerun to get cross" $*.log ) do true ; \
	done

%.bbl : %.tex biblio.bib
	pdflatex $*.tex
	-bibtex $*
	
clean: 
	bash -c "rm -f manual{.log,.bbl,.idx,.aux,.blg,Notes.bib}"
