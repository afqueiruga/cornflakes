SRCDIR = .
BUILDDIR = build
INSTALLDIR = ..

ARCH = linuxbox
ifeq ($(ARCH),linuxbox)
CC = gcc
SWIG = swig3.0
CFLAGS = -c -fpic #-framework Python
IFLAGS = -I/usr/include/python2.7/ -I/usr/include/sys/ -I$(SRCDIR)
LDSETTING = -shared
else
CC = gcc-mp-5
SWIG = swig
IFLAGS = -I/opt/local/Library/Frameworks/Python.framework/Headers/ -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/
LFLAGS = -lpython -lm
LDSETTING = -dynamiclib
endif



SWIGINFS = cornflakes_library_swig.i
WRAPPERS = $(SWIGINFS:.i=_wrap.c)
#SWIGPYS = $(SWIGINFS:.i=.py)
SRCS = hypergraph.c spatialhash.c graphers.c assemble.c kernel.c dofmap.c sample_peri.c
HEADERS = $(SRCS:.c=.h)
OBJS = $(SRCS:.c=.o)
WRAPOBJS = $(WRAPPERS:.c=.o)
TARGETS = _cornflakes_library.so 
PYTARGETS = cornflakes_library.py

TARGETSI = $(TARGETS:%=$(INSTALLDIR)/%)
PYTARGETSI = $(PYTARGETS:%=$(INSTALLDIR)/%)
WRAPPERSB = $(WRAPPERS:%=$(BUILDDIR)/%)
WRAPOBJSB = $(WRAPOBJS:%=$(BUILDDIR)/%)
OBJSB = $(OBJS:%=$(BUILDDIR)/%)

all: $(TARGETSI) $(PYTARGETSI)

$(TARGETSI): $(OBJSB) $(WRAPOBJSB)
	$(CC) $(LDSETTING) $(OBJSB) $(WRAPOBJSB) $(LFLAGS) -o $@

$(WRAPPERSB): $(BUILDDIR)/%_wrap.c : %.i $(HEADERS)
	$(SWIG) -python -outdir $(INSTALLDIR) -o $@ $<

$(OBJSB): $(BUILDDIR)/%.o : %.c %.h
	$(CC) $(CFLAGS) $(IFLAGS) $< -o $@
$(WRAPOBJSB): %.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) $(IFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(TARGETSI) $(PYTARGETSI) $(OBJSB) $(WRAPPERSB)
