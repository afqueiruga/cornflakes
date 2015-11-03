CC = gcc
SWIG = swig3.0

CFLAGS = -c -fpic
IFLAGS = -I/usr/include/python2.7/

SWIGINFS = hypergraph_swig.i
WRAPPERS = $(SWIGINFS:.i=_wrap.c)
SWIGPYS = $(SWIGINFS:.i=.py)
SRCS = hypergraph.c SpatialHash.c Graphers.c $(WRAPPERS)
OBJS = $(SRCS:.c=.o)

TARGETS = _mylibrary.so

all: $(TARGETS)

$(TARGETS): $(OBJS)
	$(CC) -shared $(OBJS) -o $@

$(WRAPPERS): %_wrap.c : %.i
	$(SWIG) -python $<

$(OBJS): %.o : %.c
	$(CC) $(CFLAGS) $(IFLAGS) $< -o $@

clean:
	rm $(TARGETS) $(OBJS) $(SWIGPYS) $(WRAPPERS)
