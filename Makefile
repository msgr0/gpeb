.PHONY: all clean test

CC  = gcc
CFLAGS = -Wall -Wextra -std=c11 -Ilibs -Ilibs/gfatools
# -march=native only where the compiler accepts it: gcc 15 on an Apple M4 Linux
# VM (OrbStack) refuses it outright ("no support for 'sme' without 'sve2'").
# Portable release builds pass NATIVE= on the command line.
NATIVE := $(shell $(CC) -march=native -E -x c /dev/null >/dev/null 2>&1 && echo -march=native)
RELEASE_FLAGS = -O2 $(NATIVE) -DNDEBUG
LDFLAGS =

GFATOOLS_LIB = libs/gfatools/libgfa1.a
GFATOOLS_DIR = libs/gfatools

# Default target
all: gpeb

ZLIB_PREFIX = $(shell brew --prefix zlib 2>/dev/null)
ZLIB_LDFLAGS = $(if $(ZLIB_PREFIX),-L$(ZLIB_PREFIX)/lib)

# The compile line is stamped in .build-flags, rewritten (so its mtime moves) only when
# it changes: a new NATIVE=, LDFLAGS= or CC rebuilds gpeb though no source file changed.
BUILD_CMD = $(CC) $(RELEASE_FLAGS) $(CFLAGS) $(LDFLAGS) $(ZLIB_LDFLAGS)
$(shell echo '$(BUILD_CMD)' | cmp -s - .build-flags 2>/dev/null || echo '$(BUILD_CMD)' > .build-flags)

# Standalone single-file guardian pebbler
gpeb: gpeb.c $(GFATOOLS_LIB) .build-flags
	$(CC) $(RELEASE_FLAGS) $(CFLAGS) $(LDFLAGS) -o gpeb gpeb.c $(GFATOOLS_LIB) $(ZLIB_LDFLAGS) -lz -lpthread -lm

# gpeb's self-test: the GPEB_TEST build is a separate binary that runs and exits.
test: gpeb.c $(GFATOOLS_LIB)
	$(CC) -O2 -DGPEB_TEST $(CFLAGS) $(LDFLAGS) -o gpeb-test gpeb.c $(GFATOOLS_LIB) $(ZLIB_LDFLAGS) -lz -lpthread -lm
	./gpeb-test

# Build gfatools library (rebuilt when any of its sources changes)
$(GFATOOLS_LIB): $(wildcard $(GFATOOLS_DIR)/*.c $(GFATOOLS_DIR)/*.h $(GFATOOLS_DIR)/gfatools/*.h)
	$(MAKE) -C $(GFATOOLS_DIR) libgfa1.a CPPFLAGS="-ffile-prefix-map=$(CURDIR)=. $(if $(ZLIB_PREFIX),-I$(ZLIB_PREFIX)/include)"

# Clean build artifacts
clean:
	rm -f gpeb gpeb-test .build-flags
	$(MAKE) -C $(GFATOOLS_DIR) clean
