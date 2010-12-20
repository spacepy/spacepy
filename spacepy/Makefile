CC=gcc
CFLAGS=-c -Wall
LDFLAGS=

OBJECTS=poppy.o randomkit.o
HEADERS=randomkit.h
LIB_UNVER=libpoppy.so
LIBRARY=libpoppy.so.0

ifdef DEBUG
CFLAGS+=-g -DDEBUG -O0
LDFLAGS+=-g -DDEBUG
else
CFLAGS+=-O2
endif

ifdef THREADED
CFLAGS+=-pthread -DTHREADED
LDFLAGS+=-pthread
endif

CFLAGS+=-fPIC
LDFLAGS+=-shared -Wl,-soname,libpoppy.so.0 -lc -lm

.PHONY: all clean

all: $(LIB_UNVER)

clean:
	rm -f *~ $(OBJECTS) $(LIBRARY) $(LIB_UNVER)

$(OBJECTS): %.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) $< -o $@

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(LIB_UNVER): $(LIBRARY)
	rm -f $(LIB_UNVER)
	ln -s $(LIBRARY) $(LIB_UNVER)
