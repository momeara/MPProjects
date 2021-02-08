
all: notebook.html notebook.pdf

R_OPTS=--no-save --no-restore --no-init-file --no-site-file # vanilla, but with --environ

notebook.html: notebook.Rmd
	R ${R_OPTS} -e "rmarkdown::render('notebok.Rmd')"

notebook.pdf: knitr_example.Rnw
	R ${R_OPTS} -e "library(knitr);knit('notebook.Rnw')"
	pdflatex notebook.tex

clean:
	rm knitr_example.tex


