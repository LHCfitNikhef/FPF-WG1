default:
	pdflatex FPF-WG1
#	bibtex FPF-WG1
	pdflatex FPF-WG1	

clean:
	rm -rf *~ *.out *.toc *.log FPF-WG1.pdf *.aux *.blg  scripts/*.pdf scripts/*~ *.tgz tables/*~

arxiv:
	rm -rf *.tgz
	tar cvzf FPF-WG1.tgz *.tex plots FPF-WG1.bbl tables
	mv FPF-WG1.tgz /Users/juanrojo/Desktop
