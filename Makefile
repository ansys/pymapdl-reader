# Simple makefile to simplify repetitive build env management tasks under posix

CODESPELL_DIRS ?= ./
CODESPELL_SKIP ?= "*.pyc,*.txt,*.gif,*.png,*.jpg,*.js,*.html,*.doctree,*.ttf,*.woff,*.woff2,*.eot,*.mp4,*.inv,*.pickle,*.ipynb,flycheck*,./.git/*,./.hypothesis/*,*.yml,./doc/build/*,./doc/images/*,./dist/*,*~,.hypothesis*,./doc/source/examples/*,*cover,*.dat,*.mac,\#*,build,./docker/mapdl/v211,./factory/*,./ansys/mapdl/reader/cython/_reader.c,./ansys/mapdl/reader/cython/_binary_reader.cpp,./ansys/mapdl/reader/cython/_cellqual.c,,./ansys/mapdl/reader/cython/_relaxmidside.c,*.inp,_*.c*,*.mypy_cache/*,*.cdb,./doc/source/examples/*"
CODESPELL_IGNORE ?= "ignore_words.txt"

all: doctest

doctest: codespell

codespell:
	@echo "Running codespell"
	@codespell $(CODESPELL_DIRS) -S $(CODESPELL_SKIP) -I $(CODESPELL_IGNORE)

pydocstyle:
	@echo "Running pydocstyle"
	@pydocstyle ansys.mapdl.reader

doctest-modules:
	@echo "Running module doctesting"
	pytest -v --doctest-modules ansys.mapdl.reader

coverage:
	@echo "Running coverage"
	@pytest -v --cov ansys.mapdl.reader

coverage-xml:
	@echo "Reporting XML coverage"
	@pytest -v --cov ansys.mapdl.reader --cov-report xml

coverage-html:
	@echo "Reporting HTML coverage"
	@pytest -v --cov ansys.mapdl.reader --cov-report html
