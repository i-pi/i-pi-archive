.PHONY : all clean
all : create_man help_files/manual.pdf

create_man:
	python create_man.py

help_files/manual.pdf: 
	lyx --export pdf2 manual.lyx
	mv manual.pdf help_files/manual.pdf
	
clean: 
	bash -c "rm -rf help_files auxiliary"
