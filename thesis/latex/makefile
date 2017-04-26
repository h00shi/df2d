bib:
	bibtex main
main:
	xelatex main.tex
clean:
	rm -f *.aux *.lof *.out *.toc *.lot *~ *bbl *log *pdf *blg *brf *dvi
all:
	xelatex main.tex
	bibtex main
	xelatex main.tex
	xelatex main.tex
