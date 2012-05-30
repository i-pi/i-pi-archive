.PHONY : all clean
all : create_man manual.pdf

create_man :
	python create_man.py

manual.pdf : 
	lyx --export pdf2 manual.lyx

%.bbl : %.tex biblio.bib
	pdflatex $*.tex
	-bibtex $*
	
clean: 
	bash -c "rm -f *.tex help.xml manual.pdf"
