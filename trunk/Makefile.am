SUBDIRS = src bindings examples

EXTRA_DIST = $(subdirs)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pHash.pc

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: @LIBTOOL_DEPS@
	$(SHELL) ./config.status --recheck
