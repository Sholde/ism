# Compilation
CC=gcc
CFLAGS=-Wall -Wextra
OFLAGS=-O3 -march=native -mtune=native # -Ofast -funroll-loops -finline-functions -ftree-vectorize
DFLAGS=-g
LFLAGS=-lm
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
	$(Q) $(LINKER) $(LFLAGS) $^ -o $@
	@if [ "$(Q)" == "@" ] ; then \
		echo "Linking complete!" ; \
		echo "Creating a binary in "$@ ; \
	fi

$(OBJDIR)/main.o: $(SRCDIR)/main.c $(VELOCITY_VERLET) $(LENNARD_JONES) $(COMMON) $(HELPER)
	$(Q) $(CC) -c $(CFLAGS) $(OFLAGS) $(DFLAGS) $(WFLAGS) $< -o $@
	@if [ "$(Q)" == "@" ] ; then \
		echo "Compiled "$<" successfully!" ; \
	fi

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(SRCDIR)/%.h
	$(Q) $(CC) -c $(CFLAGS) $(OFLAGS) $(DFLAGS) $(WFLAGS) $< -o $@
	@if [ "$(Q)" == "@" ] ; then \
		echo "Compiled "$<" successfully!" ; \
	fi

# Dependencies variable
VELOCITY_VERLET= $(SRCDIR)/velocity_verlet.c $(SRCDIR)/velocity_verlet.h
LENNARD_JONES= $(SRCDIR)/lennard_jones.c $(SRCDIR)/lennard_jones.h
COMMON= $(SRCDIR)/common.c $(SRCDIR)/common.h
HELPER= $(SRCDIR)/helper.h

# Dependencies target
$(SRCDIR)/velocity_verlet.c: $(LENNARD_JONES) $(COMMON) $(HELPER)

$(SRCDIR)/lennard_jones.c: $(COMMON) $(HELPER)

$(SRCDIR)/common.c: $(HELPER)

# Cleanup
clean:
	$(Q) rm -Rf *~ **/*~ $(OBJDIR) $(BINDIR)
	@if [ "$(Q)" == "@" ] ; then \
		echo "Cleanup complete!" ; \
	fi
