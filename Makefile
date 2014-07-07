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
