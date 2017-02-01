LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
bin_PROGRAMS = otc-tol-ws


LIBS = $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS) -lrestbed -lssl -lcurl -lcrypto

otc_tol_ws_SOURCES = tolws.cpp
otc_tol_ws_CPPFLAGS = $(AM_CPPFLAGS)

