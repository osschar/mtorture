AUTO_FILES := $(addprefix Matrix, $(addsuffix .ah, 33x33 66x66 Sym33x33 Sym66x66))

all: auto

auto: ${AUTO_FILES}

${AUTO_FILES}:
	./gen-matmul.pl

clean:
	rm -f *.ah

distclean: clean

echo:
	@echo AUTO_FILES = ${AUTO_FILES}
