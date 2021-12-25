# Compilation
CC=gcc
CFLAGS=-Wall -Wextra
OFLAGS=-O3 -march=native -mtune=native
DFLAGS=-g
LFLAGS=
WFLAGS=-Wno-incompatible-pointer-types

# Linking
TARGET=main
LINKER=$(CC)

# Diretories
SRCDIR=src
OBJDIR=obj
BINDIR=bin

# Get filenames
SRCS=$(wildcard $(SRCDIR)/*.c)
OBJS=$(SRCS:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

# General
Q=@

# Phony
.PHONY: all clean

# Target
all: dir $(BINDIR)/$(TARGET)

dir:
	@mkdir -p src obj bin

$(BINDIR)/$(TARGET): $(OBJS)
	$(Q) $(LINKER) $(LFLAGS)$^ -o $@
	@if [ "$(Q)" == "@" ] ; then \
		echo "Linking complete!" ; \
		echo "Creating a binary in "$@ ; \
	fi

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(Q) $(CC) -c $(CFLAGS) $(OFLAGS) $(DFLAGS) $(WFLAGS) $< -o $@
	@if [ "$(Q)" == "@" ] ; then \
		echo "Compiled "$<" successfully!" ; \
	fi

$(SRCDIR)/main.c: $(SRCDIR)/lennard_jones.c $(SRCDIR)/lennard_jones.h $(SRCDIR)/common.c $(SRCDIR)/common.h $(SRCDIR)/helper.h

$(SRCDIR)/lennard_jones.c: $(SRCDIR)/lennard_jones.h $(SRCDIR)/common.c $(SRCDIR)/common.h $(SRCDIR)/helper.h

$(SRCDIR)/common.c: $(SRCDIR)/common.h $(SRCDIR)/helper.h

clean:
	$(Q) rm -Rf *~ **/*~ $(OBJDIR) $(BINDIR)
	@if [ "$(Q)" == "@" ] ; then \
		echo "Cleanup complete!" ; \
	fi
