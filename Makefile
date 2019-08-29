CXX=clang++

.NOTPARALLEL:
all: run clean build-doc 

build-doc:
	pandoc --filter pandoc-crossref --from=markdown --to=latex --standalone --variable=colorlinks:true -V fontsize=12pt readme.md -o readme.pdf

run: build
	./test

build:
	${CXX} -std=c++14 test.cpp -I/usr/local/include/eigen3 -I/usr/include/eigen3 -o test

clean:
	rm -f test readme.pdf
