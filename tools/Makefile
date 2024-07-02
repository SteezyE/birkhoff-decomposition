CC=g++
CCFLAGS=-march=native -Ofast
TARGETS=mm_generator_uniform mm_generator_skewed

all: $(TARGETS)

mm_generator_uniform: mm_generator_uniform.cpp
	$(CC) $(CCFLAGS) -o $@ $^

mm_generator_skewed: mm_generator_skewed.cpp
	$(CC) $(CCFLAGS) -o $@ $^

.PHONY: clean
clean: 
	rm -rf $(TARGETS)
