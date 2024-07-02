CC=gcc
CCFLAGS=-march=native -Ofast -lm -flto -Wno-alloc-size-larger-than
TARGETS=birkhoff beta_solstice binary_bottleneck kpack birkhoff.o beta_solstice.o \
        binary_bottleneck.o kpack.o bottleneck_matching.o uno.o mmio.o cscio.o cheap.o matching.o

all: $(TARGETS)

birkhoff: birkhoff.o mmio.o cheap.o matching.o cscio.o
	$(CC) $(CCFLAGS) -o $@ $^

binary_bottleneck: binary_bottleneck.o bottleneck_matching.o mmio.o cheap.o matching.o cscio.o
	$(CC) $(CCFLAGS) -o $@ $^

beta_solstice: beta_solstice.o mmio.o cheap.o matching.o cscio.o
	$(CC) $(CCFLAGS) -o $@ $^

kpack: kpack.o bottleneck_matching.o uno.o mmio.o cheap.o matching.o cscio.o
	$(CC) $(CCFLAGS) -o $@ $^

birkhoff.o: birkhoff.c
	$(CC) $(CCFLAGS) -c $<

binary_bottleneck.o: binary_bottleneck.c
	$(CC) $(CCFLAGS) -c $<

beta_solstice.o: beta_solstice.c
	$(CC) $(CCFLAGS) -c $<

kpack.o: kpack.c
	$(CC) $(CCFLAGS) -c $<

bottleneck_matching.o: bottleneck_matching.c
	$(CC) $(CCFLAGS) -c $<

uno.o: extern/perfect-matching-enumeration/uno.c
	$(CC) $(CCFLAGS) -c $<

mmio.o: extern/mmio/mmio.c extern/mmio/mmio.h
	$(CC) $(CCFLAGS) -c $<

cscio.o: cscio.c cscio.h
	$(CC) $(CCFLAGS) -c $<

cheap.o: extern/matchmaker/cheap.c extern/matchmaker/matchmaker.h
	$(CC) $(CCFLAGS) -c $<

matching.o: extern/matchmaker/matching.c extern/matchmaker/matchmaker.h
	$(CC) $(CCFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(TARGETS)
