# MIT License
#
# Copyright (c) 2022 André Bianchessi
# Copyright (c) 2019-2022 Carsten Gips
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

###############################################################################
## Setup (public)
###############################################################################


META = metadata.yaml

## Add chapters here inh the order they appear with
## SRC += chapters/XX_chapter.md
SRC  = chapters/acronyms.md
SRC += chapters/01_introduction.md
SRC += chapters/02_objectives.md
SRC += chapters/03_literature_review.md
SRC += chapters/04_methods.md
SRC += chapters/05_results.md
SRC += chapters/06_discussion.md

BIBFILE = references.bib
CSL_FILE = ieee.csl

APPENDIX = chapters/appendix.md


###############################################################################
## Internal setup (do not change)
###############################################################################


## Auxiliary files
## (Do not change!)
TITLEPAGE   = chapters/titlepage.tex
FRONTMATTER = chapters/frontmatter.tex
BACKMATTER  = chapters/backmatter.tex
REFERENCES  = chapters/references.md

TMP1 = $(TITLEPAGE:%.tex=%.filled.tex)
TMP2 = $(FRONTMATTER:%.tex=%.filled.tex)
TMP3 = $(BACKMATTER:%.tex=%.filled.tex)
TMP  = $(TMP1) $(TMP2) $(TMP3)


## Pandoc options
AUX_OPTS  = --wrap=preserve
AUX_OPTS += -M cleanthesis=true -M cleanthesisbibfile=$(BIBFILE:%.bib=%)

OPTIONS  = -f markdown
OPTIONS += --pdf-engine=pdflatex
OPTIONS += --standalone

OPTIONS += --metadata-file=$(META)

OPTIONS += --include-in-header=$(TMP1)
OPTIONS += --include-before-body=$(TMP2)
OPTIONS += --include-after-body=$(TMP3)

OPTIONS += --filter pandoc-crossref
OPTIONS += --citeproc
OPTIONS += -M bibliography=$(BIBFILE)
OPTIONS += -M link-citations=true

OPTIONS += --listings

OPTIONS += -V documentclass=scrbook
OPTIONS += -V papersize=a4
OPTIONS += -V fontsize=11pt

OPTIONS += -V classoption:open=right
OPTIONS += -V classoption:twoside=true
OPTIONS += -V classoption:cleardoublepage=empty
OPTIONS += -V classoption:clearpage=empty

OPTIONS += -V geometry:top=30mm
OPTIONS += -V geometry:left=25mm
OPTIONS += -V geometry:bottom=30mm
OPTIONS += -V geometry:width=150mm
OPTIONS += -V geometry:bindingoffset=6mm

OPTIONS += --toc
OPTIONS += --toc-depth=3
OPTIONS += --number-sections

OPTIONS += -V colorlinks=true
OPTIONS += --include-in-header=include-header.tex $(AUX_OPTS)

WEBSITE_OPTIONS  = -f markdown
WEBSITE_OPTIONS += --template template.html
WEBSITE_OPTIONS += --toc --standalone
WEBSITE_OPTIONS += --toc-depth=3
WEBSITE_OPTIONS += --metadata-file=metadata.yaml
WEBSITE_OPTIONS += --filter pandoc-crossref
WEBSITE_OPTIONS += --citeproc
WEBSITE_OPTIONS += -M bibliography=$(BIBFILE)
WEBSITE_OPTIONS += -M link-citations=true
WEBSITE_OPTIONS += --mathjax

###############################################################################
## Targets
###############################################################################
.DEFAULT_GOAL := pdf

## Build auxiliary files (title page, frontmatter, backmatter, references)
$(TMP): %.filled.tex: %.tex $(META)
	pandoc $(AUX_OPTS) --template=$< --metadata-file=$(META) -o $@ $<

## Build thesis
pdf: $(SRC) $(REFERENCES) $(APPENDIX) $(META) $(BIBFILE) $(TMP)
	pandoc ${OPTIONS} -o thesis.pdf $(SRC) $(REFERENCES) --csl $(CSL_FILE) $(APPENDIX)

website: $(SRC) $(REFERENCES) template.html
	pandoc $(WEBSITE_OPTIONS) -o website/index.html $(SRC) $(REFERENCES) --csl $(CSL_FILE)

## Clean-up: Remove thesis and temporary (generated) files
clean:
	rm -rf thesis.pdf website/index.html $(TMP)

