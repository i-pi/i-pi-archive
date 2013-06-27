.PHONY : all clean aux
all : manual.pdf

aux:
	python create_man.py

manual.pdf: aux
	lyx --export pdf2 manual.lyx
	mv manual.pdf help_files	
	
clean: 
	bash -c "rm -rf help_files input_docs"
