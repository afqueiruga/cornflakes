ARCH = linuxbox
ifeq ($(ARCH),linuxbox)
CC = gcc
SWIG = swig3.0
CFLAGS = -c -fpic #-framework Python
IFLAGS = -I/usr/include/python2.7/ -I/usr/include/sys/
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
SWIGPYS = $(SWIGINFS:.i=.py)
SRCS = hypergraph.c spatialhash.c graphers.c assemble.c kernel.c $(WRAPPERS)
OBJS = $(SRCS:.c=.o)

TARGETS = _cornflakes_library.so

all: $(TARGETS)

$(TARGETS): $(OBJS)
	$(CC) $(LDSETTING) $(OBJS) $(LFLAGS) -o $@

$(WRAPPERS): %_wrap.c : %.i
	$(SWIG) -python $<

$(OBJS): %.o : %.c
	$(CC) $(CFLAGS) $(IFLAGS) $< -o $@

clean:
	rm $(TARGETS) $(OBJS) $(SWIGPYS) $(WRAPPERS)
