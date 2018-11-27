.PHONY: all clean realclean ALWAYS

TARGETS:=laplacian-spectrum.pdf

# Let latexmk handle tex dependencies,
# but also generate a dependency file for generated images etc.

LATEXMK:=latexmk
LATEXMKOPT:=-e '$$clean_ext=q/run.xml/' -bibtex -xelatex -use-make -interaction=nonstopmode

ifdef force
    LATEXMKOPT+=-gg
endif

ifndef errstop
    LATEXMKOPT+=-f
endif

all: $(TARGETS)

%.pdf: %.tex ALWAYS
	$(LATEXMK) $(LATEXMKOPT) $< --deps-out="$@.d"

$(foreach t,$(TARGETS),$(eval -include $(t).d))

clean:
	for t in $(TARGETS); do $(LATEXMK) $(LATEXMKOPT) -c "$$t"; done

realclean:
	for t in $(TARGETS); do \
	$(LATEXMK) $(LATEXMKOPT) -C "$$t"; \
	rm -f "$$t.d"; \
	done
