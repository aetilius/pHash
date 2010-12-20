include_HEADERS = pHash-config.h
SUBDIRS = src bindings examples
EXTRA_DIST = $(subdirs) libpHash.spec

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pHash.pc

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: @LIBTOOL_DEPS@
	$(SHELL) ./config.status --recheck

dist-hook:
	rm -rf `find $(distdir) -name \.svn`
