clean:
	rm -f -r build/
	rm -f src/python/cython_utils.so

clean-cpp:
	rm -f src/python/cython_utils.cpp

clean-all: clean clean-cpp

.PHONY: build
build: clean
	python setup.py build_ext --inplace

.PHONY: cython-build
cython-build: clean clean-cpp
	python setup.py build_ext --inplace --use-cython

.PHONY: tests
tests:
	hash nosetests 2>/dev/null || { echo -e >&2 "############################\nI require the python library \"nose\" for unit tests but it's not installed.  Aborting.\n############################\n"; exit 1; } 
	nosetests --logging-level=INFO tests/
