# ----------------------- shared settings

SHARED_CFLAGS = -g -DPD -g -DPD -Wall -W -Wshadow -Wstrict-prototypes \
    -Wno-unused -Wno-parentheses -Wno-switch

# ----------------------- LINUX i386 -----------------------

pd_linux: oscbank~.pd_linux

.SUFFIXES: .pd_linux

LINUXCFLAGS = -DPD -O2 -funroll-loops -fomit-frame-pointer \
    -Wall -W -Wshadow -Wstrict-prototypes -Werror \
    -Wno-unused -Wno-parentheses -Wno-switch

LINUXINCLUDE =  -I /usr/local/src/pd

.c.pd_linux:
	cc $(LINUXCFLAGS) $(LINUXINCLUDE) -o $*.o -c $*.c
	ld --export-dynamic  -shared -o $*.pd_linux $*.o -lc -lm 
	strip --strip-unneeded $*.pd_linux
	rm $*.o


# ----------------------- Mac OSX -----------------------

pd_darwin: oscbank~.pd_darwin

.SUFFIXES: .pd_darwin

DARWINCFLAGS = #(SHARED_CFLAGS)

DARWIN_INCLUDE = -I /Applications/Pd-extended.app/Contents/Resources/include/

.c.pd_darwin:
	cc -arch i386 $(DARWINCFLAGS) $(DARWIN_INCLUDE) -o $*.o -c $*.c
	cc -arch i386 -bundle -undefined suppress -flat_namespace -o $*.pd_darwin $*.o
	rm -f $*.o

clean:
	rm *.pd_darwin
# ----------------------- NT -----------------------

pd_nt: oscbank~.dll

.SUFFIXES: .obj .dll

PDNTCFLAGS = /W3 /WX /DNT /DPD /nologo
VC="D:\Program Files\Microsoft Visual Studio\Vc98"

PDNTINCLUDE = /I. /I\tcl\include /I..\..\src /I$(VC)\include

PDNTLDIR = $(VC)\lib
PDNTLIB = $(PDNTLDIR)\libc.lib \
	$(PDNTLDIR)\oldnames.lib \
	$(PDNTLDIR)\kernel32.lib \
	..\..\bin\pd.lib 

.c.dll:
	cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	link /dll /export:$*_setup $*.obj $(PDNTLIB)

# override explicitly for tilde objects like this:
dspobj~.dll: dspobj~.c; 
	cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	link /dll /export:dspobj_tilde_setup $*.obj $(PDNTLIB)