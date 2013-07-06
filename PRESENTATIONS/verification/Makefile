TEX_TARGETS = 
TARGETS = output presentation.pdf $(TEX_TARGETS)

presentation.pdf: $(TEX_TARGETS)

%.pdf: %.tex bibliography.bib
	pdflatex --enable-write18 $<
	if fgrep "Rerun to" $*.log; then $(MAKE) --assume-new $^ $@;fi
ifndef NOBIBTEX
	if fgrep "There were undefined"  $*.log;\
	then bibtex $*; \
	$(MAKE) NOBIBTEX=true --assume-new $^ $@;fi
endif

.PHONY: clean

# We don't rm *.ps because that might clobber image files. Instead, put any
# ps files which should go in $(PSFILES)
clean: 
	rm -r $(TARGETS) $(PSFILES) *.dvi *.aux *.log *.bbl *.blg \
*.out *.snm *.toc *.nav *.vrb 2>/dev/null||true
