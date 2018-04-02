SOURCES := main.cpp DevelopingCell.cpp Cell.cpp
# Objs are all the sources, with .cpp replaced by .o
OBJS := $(SOURCES:.cpp=.o)

all: main

# Compile the binary 't' by calling the compiler with cflags, lflags, and any libs (if defined) and the list of objects.
main: $(OBJS)
	$(CC) $(CFLAGS) -o t $(OBJS) $(LFLAGS) $(LIBS)

# Get a .o from a .cpp by calling compiler with cflags and includes (if defined)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

