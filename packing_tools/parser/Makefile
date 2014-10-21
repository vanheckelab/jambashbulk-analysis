ifeq ($(OS),Windows_NT)
    dllname = _parser.dll
    CCFLAGS += -D WIN32
else
    dllname = _parser.so
    CCFLAGS += -D LINUX
endif
dll = $(dllname)
srcfiles = parser.c

all:
	gcc $(CCFLAGS) -fPIC -shared -Wl,-soname,$(dllname) -o $(dll) $(srcfiles)
