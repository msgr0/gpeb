# Standalone .pex tools, one C file each, no libraries:  make   /  make clean
CC     ?= cc
CFLAGS ?= -O2 -Wall
BINS    = pex_rare pex_bycolour

all: $(BINS)

# make's built-in %: %.c rule links each binary from its own source
clean:
	rm -f $(BINS)

.PHONY: all clean
