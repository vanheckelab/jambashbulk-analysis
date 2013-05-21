dllname = _parser.so
dll = $(dllname)
srcfiles = parser.c

all:
	gcc -fPIC -shared -Wl,-soname,$(dllname) -o $(dll) $(srcfiles)
