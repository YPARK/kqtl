RMD := $(wildcard *.rmd)
HTML := $(RMD:.rmd=.html)
PDF := $(RMD:.rmd=.pdf)
DOCX := $(RMD:.rmd=.docx)

all: $(HTML)
docx: $(DOCX)
pdf: $(PDF)

%.html: %.rmd
	Rscript -e "library(rmarkdown); render('$<');"

%.pdf: %.rmd
	Rscript -e "library(rmarkdown); render('$<', 'pdf_document');"

%.docx: %.rmd
	Rscript -e "library(rmarkdown); render('$<', 'word_document');"

