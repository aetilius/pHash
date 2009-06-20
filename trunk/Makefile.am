SUBDIRS = src bindings examples

DIST_SUBDIRS = $(subdirs)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pHash.pc

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: @LIBTOOL_DEPS@
	$(SHELL) ./config.status --recheck
