# GUINEA-PIG Makefile

prefix = /home/baonguye/guinea-pig/install
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
libdir = ${exec_prefix}/lib

all: guinea

guinea_nofftw:
	$(MAKE) -C src guinea_nofftw 2>&1 | tee log.guinea-nofftw.$$$$

guinea:
	$(MAKE) -C src guinea 2>&1 | tee log.guinea.$$$$

guinearoot: 
	$(MAKE) -C src guinearoot 2>&1 | tee log.guinea-root.$$$$

install: all
	test -d $(bindir) || mkdir -p $(bindir)
	cd src && cp guinea $(bindir) && cd ..
clean: 
	rm -f log.*
	$(MAKE) -C src clean

distclean: clean
	rm -f config.log config.status Makefile src/Makefile src/config.h src/stamp-h
