TEX_TARGETS = branching.tex branch_merging.tex dvcs.tex vcs.tex \
pull_request_2.tex pull_request.tex
TARGETS = output presentation.pdf $(TEX_TARGETS)

presentation.pdf: $(TEX_TARGETS)

$(TEX_TARGETS): vcs.py
	python vcs.py $(basename $@)

%.pdf: %.tex 
	pdflatex --enable-write18 $<
	if fgrep "Rerun to" $*.log; then $(MAKE) --assume-new $^ $@;fi

.PHONY: clean

# We don't rm *.ps because that might clobber image files. Instead, put any
# ps files which should go in $(PSFILES)
clean: 
	rm -r $(TARGETS) $(PSFILES) *.dvi *.aux *.log *.bbl *.blg \
*.out *.snm *.toc *.nav *.vrb 2>/dev/null||true
