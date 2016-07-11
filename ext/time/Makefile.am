bin_PROGRAMS = time
info_TEXINFOS = time.texi
time_SOURCES = time.c getopt.c getopt1.c error.c resuse.c version.c \
	getopt.h getpagesize.h resuse.h port.h wait.h
EXTRA_DIST = texinfo.tex

DISTCLEANFILES = stamp-v version.c
version.c: stamp-v
stamp-v: Makefile
	rm -f t-version.c
	echo 'const char *version_string = "'GNU @PACKAGE@ @VERSION@'";' \
		> t-version.c
	if cmp -s version.c t-version.c; then	\
	  rm t-version.c;			\
	else					\
	  mv t-version.c version.c;		\
	fi
	echo timestamp > $@
