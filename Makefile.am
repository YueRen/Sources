ACLOCAL_AMFLAGS = -I m4

SUBDIRS = resources omalloc factory \
libpolys gfanlib IntegerProgramming \
kernel Singular \
dox emacs redhat desktop

EXTRA_DIST = README README.md README.pkg autogen.sh doxy redhat/singular.spec.in GPL2 GPL3

bin_SCRIPTS = libsingular-config

install-data-local:
	if test -e $(srcdir)/doc/doc.tbz2; then\
	  (cat ${top_srcdir}/doc/doc.tbz2| (cd $(DESTDIR)$(datadir)/singular; tar jxf -));\
	   mkdir $(DESTDIR)$(datadir)/info;\
	   mkdir $(DESTDIR)$(datadir)/doc;\
	   mv  $(DESTDIR)$(datadir)/singular/singular.hlp  $(DESTDIR)$(datadir)/info/.;\
	fi

uninstall-local:
	-(cd $(DESTDIR)$(datadir); rm -rf singular info/singular.hlp)

configheaderdir = ${includedir}/singular
nodist_configheader_HEADERS = singularconfig.h
DISTCLEANFILES = $(nodist_configheader_HEADERS)

dist_man1_MANS = doc/Singular.man doc/ESingular.man doc/TSingular.man

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = Singular.pc
