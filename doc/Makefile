# Makefile for Sphinx documentation
#

PYVER = 3.5
PYTHON = python$(PYVER)

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = LANG=C sphinx-build
PAPER         =

FILES=

#BUILD_DIR 		= /mnt/data1/build/sphinx-build/scikit-nano
BUILD_DIR 		= build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILD_DIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck \
        dist dist-build

#------------------------------------------------------------------------------

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  html-sknanowiki  to make standalone HTML files with scipy.org theming"
	@echo "  pickle    to make pickle files (usable by e.g. sphinx-web)"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  dist PYVER=... to make a distribution-ready tree"
	@echo "  upload USER=...  to upload results to docs.scipy.org"

clean:
	-rm -rf $(BUILD_DIR)/* source/generated *.xyz


#------------------------------------------------------------------------------
# Automated generation of all documents
#------------------------------------------------------------------------------

# Build the current scipy version, and extract docs from it.
# We have to be careful of some issues:
#
# - Everything must be done using the same Python version
# - We must use eggs (otherwise they might override PYTHONPATH on import).
# - Different versions of easy_install install to different directories (!)
#

INSTALL_DIR = $(CURDIR)/$(BUILD_DIR)/inst-dist/
INSTALL_PPH = $(INSTALL_DIR)/lib/python$(PYVER)/site-packages:$(INSTALL_DIR)/local/lib/python$(PYVER)/site-packages:$(INSTALL_DIR)/lib/python$(PYVER)/dist-packages:$(INSTALL_DIR)/local/lib/python$(PYVER)/dist-packages

DIST_VARS=PYTHON="PYTHONPATH=$(INSTALL_PPH):$$PYTHONPATH python$(PYVER)" SPHINXBUILD="LANG=C PYTHONPATH=$(INSTALL_PPH):$$PYTHONPATH python$(PYVER) `which sphinx-build`"

dist:
	make $(DIST_VARS) real-dist

real-dist: dist-build html html-sknanowiki
	test -d $(BUILD_DIR)/latex || make latex
	make -C $(BUILD_DIR)/latex all-pdf
	-test -d $(BUILD_DIR)/htmlhelp || make htmlhelp-build
	-rm -rf $(BUILD_DIR)/dist
	mkdir -p $(BUILD_DIR)/dist
	cp -r $(BUILD_DIR)/html-sknanowiki $(BUILD_DIR)/dist/reference
	touch $(BUILD_DIR)/dist/index.html
	(cd $(BUILD_DIR)/html && zip -9qr ../dist/scipy-html.zip .)
	cp $(BUILD_DIR)/latex/scipy*.pdf $(BUILD_DIR)/dist
	-zip $(BUILD_DIR)/dist/scipy-chm.zip $(BUILD_DIR)/htmlhelp/scipy.chm
	chmod ug=rwX,o=rX -R $(BUILD_DIR)/dist
	find $(BUILD_DIR)/dist -type d -print0 | xargs -0r chmod g+s
	cd $(BUILD_DIR)/dist && tar czf ../dist.tar.gz *

dist-build:
	rm -f ../dist/*.egg
	cd .. && $(PYTHON) setupegg.py bdist_egg
	install -d $(subst :, ,$(INSTALL_PPH))
	$(PYTHON) `which easy_install` --prefix=$(INSTALL_DIR) ../dist/*.egg


#------------------------------------------------------------------------------
# Basic Sphinx generation rules for different formats
#------------------------------------------------------------------------------

html:
	mkdir -p $(BUILD_DIR)/html $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILD_DIR)/html $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILD_DIR)/html."

html-sknanowiki:
	mkdir -p $(BUILD_DIR)/html $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -t sknanowiki -b html $(ALLSPHINXOPTS) $(BUILD_DIR)/html-sknanowiki $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILD_DIR)/html."

pickle:
	mkdir -p $(BUILD_DIR)/pickle $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) $(BUILD_DIR)/pickle $(FILES)
	@echo
	@echo "Build finished; now you can process the pickle files or run"
	@echo "  sphinx-web $(BUILD_DIR)/pickle"
	@echo "to start the sphinx-web server."

web: pickle

htmlhelp:
	mkdir -p $(BUILD_DIR)/htmlhelp $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) $(BUILD_DIR)/htmlhelp $(FILES)
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILD_DIR)/htmlhelp."

htmlhelp-build: htmlhelp $(BUILD_DIR)/htmlhelp/scipy.chm
%.chm: %.hhp
	-hhc.exe $^

latex:
	mkdir -p $(BUILD_DIR)/latex $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILD_DIR)/latex $(FILES)
	$(PYTHON) postprocess.py tex $(BUILD_DIR)/latex/*.tex
	perl -pi -e 's/\t(latex.*|pdflatex) (.*)/\t-$$1 -interaction batchmode $$2/' $(BUILD_DIR)/latex/Makefile
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILD_DIR)/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

coverage: build
	mkdir -p $(BUILD_DIR)/coverage $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b coverage $(ALLSPHINXOPTS) $(BUILD_DIR)/coverage $(FILES)
	@echo "Coverage finished; see c.txt and python.txt in $(BUILD_DIR)/coverage"

changes:
	mkdir -p $(BUILD_DIR)/changes $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) $(BUILD_DIR)/changes $(FILES)
	@echo
	@echo "The overview file is in build/changes."

linkcheck:
	mkdir -p $(BUILD_DIR)/linkcheck $(BUILD_DIR)/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) $(BUILD_DIR)/linkcheck $(FILES)
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."
