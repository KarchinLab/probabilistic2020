SRC_DIR = permutation2020

clean:
	rm -f -r build/
	rm -f $(SRC_DIR)/cython/cutils.so
	rm -f $(SRC_DIR)/cython/gaussian_kde.so
	rm -f $(SRC_DIR)/cython/uniform_kde.so
	find . -name '*.pyc' -delete

clean-cpp:
	rm -f $(SRC_DIR)/cython/cutils.cpp
	rm -f $(SRC_DIR)/cython/gaussian_kde.cpp
	rm -f $(SRC_DIR)/cython/uniform_kde.cpp

clean-all: clean clean-cpp

.PHONY: build
build: clean
	python setup.py build_ext --inplace

.PHONY: cython-build
cython-build: clean clean-cpp
	python setup.py build_ext --inplace --use-cython

.PHONY: tests
tests:
	@hash nosetests 2>/dev/null || { echo -e >&2 "############################\nI require the python library \"nose\" for unit tests but it's not installed.  Aborting.\n############################\n"; exit 1; } 
	nosetests --logging-level=INFO tests/
