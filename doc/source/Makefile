# Makefile for Sphinx documentation
#

PYVER = 3.5
PYTHON = python$(PYVER)

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = LANG=C sphinx-build
PAPER         =

BUILDDIR 		= _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml web pickle json htmlhelp latex changes linkcheck \
        dist dist-build

#------------------------------------------------------------------------------

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  html-sknanowiki  to make standalone HTML files with scikit-nano.org theming"
	@echo "  pickle    to make pickle files (usable by e.g. sphinx-web)"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  dist PYVER=... to make a distribution-ready tree"
	@echo "  upload USER=...  to upload results to docs.scikit-nano.org"

clean:
	-rm -rf $(BUILDDIR)
	-rm -rf api
	-rm *.xyz


#------------------------------------------------------------------------------
# Automated generation of all documents
#------------------------------------------------------------------------------

# Build the current scikit-nano version, and extract docs from it.
# We have to be careful of some issues:
#
# - Everything must be done using the same Python version
# - We must use eggs (otherwise they might override PYTHONPATH on import).
# - Different versions of easy_install install to different directories (!)
#

INSTALL_DIR = $(CURDIR)/$(BUILDDIR)/inst-dist/
INSTALL_PPH = $(INSTALL_DIR)/lib/python$(PYVER)/site-packages:$(INSTALL_DIR)/local/lib/python$(PYVER)/site-packages:$(INSTALL_DIR)/lib/python$(PYVER)/dist-packages:$(INSTALL_DIR)/local/lib/python$(PYVER)/dist-packages

DIST_VARS=PYTHON="PYTHONPATH=$(INSTALL_PPH):$$PYTHONPATH python$(PYVER)" SPHINXBUILD="LANG=C PYTHONPATH=$(INSTALL_PPH):$$PYTHONPATH python$(PYVER) `which sphinx-build`"

dist:
	make $(DIST_VARS) real-dist

real-dist: dist-build html html-sknanowiki
	test -d $(BUILDDIR)/latex || make latex
	make -C $(BUILDDIR)/latex all-pdf
	-test -d $(BUILDDIR)/htmlhelp || make htmlhelp-build
	-rm -rf $(BUILDDIR)/dist
	mkdir -p $(BUILDDIR)/dist
	cp -r $(BUILDDIR)/html-sknanowiki $(BUILDDIR)/dist/reference
	touch $(BUILDDIR)/dist/index.html
	(cd $(BUILDDIR)/html && zip -9qr ../dist/scikit-nano-html.zip .)
	cp $(BUILDDIR)/latex/scikit-nano*.pdf $(BUILDDIR)/dist
	-zip $(BUILDDIR)/dist/scikit-nano-chm.zip $(BUILDDIR)/htmlhelp/scikit-nano.chm
	chmod ug=rwX,o=rX -R $(BUILDDIR)/dist
	find $(BUILDDIR)/dist -type d -print0 | xargs -0r chmod g+s
	cd $(BUILDDIR)/dist && tar czf ../dist.tar.gz *

dist-build:
	rm -f ../dist/*.egg
	cd .. && $(PYTHON) setupegg.py bdist_egg
	install -d $(subst :, ,$(INSTALL_PPH))
	$(PYTHON) `which easy_install` --prefix=$(INSTALL_DIR) ../dist/*.egg


#------------------------------------------------------------------------------
# Basic Sphinx generation rules for different formats
#------------------------------------------------------------------------------

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) $(BUILDDIR)/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

html-sknanowiki:
	mkdir -p $(BUILDDIR)/html $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -t sknanowiki -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html-sknanowiki
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pickle:
	mkdir -p $(BUILDDIR)/pickle $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILDDIR)/pickle
	@echo
	@echo "Build finished; now you can process the pickle files or run"
	@echo "  sphinx-web $(BUILDDIR)/pickle"
	@echo "to start the sphinx-web server."

web: pickle

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) $(BUILDDIR)/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	mkdir -p $(BUILDDIR)/htmlhelp $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILDDIR)/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

htmlhelp-build: htmlhelp $(BUILDDIR)/htmlhelp/scikit-nano.chm
%.chm: %.hhp
	-hhc.exe $^

latex:
	mkdir -p $(BUILDDIR)/latex $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	$(PYTHON) postprocess.py tex $(BUILDDIR)/latex/*.tex
	perl -pi -e 's/\t(latex.*|pdflatex) (.*)/\t-$$1 -interaction batchmode $$2/' $(BUILDDIR)/latex/Makefile
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

coverage: build
	mkdir -p $(BUILDDIR)/coverage $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b coverage $(ALLSPHINXOPTS) $(BUILDDIR)/coverage
	@echo "Coverage finished; see c.txt and python.txt in $(BUILDDIR)/coverage"

changes:
	mkdir -p $(BUILDDIR)/changes $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILDDIR)/changes
	@echo
	@echo "The overview file is in build/changes."

linkcheck:
	mkdir -p $(BUILDDIR)/linkcheck $(BUILDDIR)/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILDDIR)/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."
