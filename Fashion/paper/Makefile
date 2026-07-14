# Build the manuscript and supplement with tectonic.
# Figures are read from ../figures/manuscript/ (committed PDFs).
TECTONIC ?= tectonic

.PHONY: all manuscript supplement cover clean

all: manuscript supplement cover

manuscript: manuscript.tex
	$(TECTONIC) manuscript.tex

supplement: supplement.tex
	$(TECTONIC) supplement.tex

cover: cover_letter_phenomenon.tex cover_letter_reproducibility.tex
	$(TECTONIC) cover_letter_phenomenon.tex
	$(TECTONIC) cover_letter_reproducibility.tex

clean:
	rm -f manuscript.pdf supplement.pdf cover_letter_phenomenon.pdf cover_letter_reproducibility.pdf *.aux *.log *.out
