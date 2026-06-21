# Build the manuscript and supplement with tectonic.
# Figures are read from ../figures/manuscript/ (committed PDFs).
TECTONIC ?= tectonic

.PHONY: all manuscript supplement cover clean

all: manuscript supplement cover

manuscript: manuscript.tex
	$(TECTONIC) manuscript.tex

supplement: supplement.tex
	$(TECTONIC) supplement.tex

cover: cover_letter.tex
	$(TECTONIC) cover_letter.tex

clean:
	rm -f manuscript.pdf supplement.pdf cover_letter.pdf *.aux *.log *.out
