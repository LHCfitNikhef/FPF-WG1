default:
	pdflatex main
#	bibtex main
	pdflatex main	

clean:
	rm -rf *~ *.out *.toc *.log main.pdf *.aux *.blg  scripts/*.pdf scripts/*~ *.tgz tables/*~

arxiv:
	rm -rf *.tgz
	tar cvzf main.tgz *.tex plots main.bbl tables
	mv main.tgz /Users/juanrojo/Desktop
