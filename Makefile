# Simple makefile to simplify repetitive build env management tasks under posix

CODESPELL_DIRS ?= ./
CODESPELL_SKIP ?= "*.pyc,*.txt,*.gif,*.png,*.jpg,*.js,*.html,*.doctree,*.ttf,*.rst,*.npz,*.eot,*.mp4,*.inv,*.pickle,*.ipynb,flycheck*,./.git/*,./.hypothesis/*,*.yml,./docs/build/*,./docs/images/*,./dist/*,*~,.hypothesis*,./docs/source/examples/*,*cover,*.dat,*.mac,\#*,build,./ansys/mapdl/reader/cython/_relaxmidside.c,*.inp,_*.c*,*.full,*.esav,*.rth,*.emat,*.npy,*.mypy_cache/*"
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
	@echo "Runnnig module doctesting"
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
