LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
bin_PROGRAMS = otc-tol-ws


LIBS = $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS) -lrestbed -lssl -lcurl -lcrypto

otc_tol_ws_SOURCES = tolws.cpp tolwsadaptors.cpp tolwsbooting.cpp nexson/nexson.cpp trees_to_serve.cpp tnrs/context.cpp tnrs/nomenclature.cpp
otc_tol_ws_CPPFLAGS = $(AM_CPPFLAGS)

check:
	python $(abs_top_srcdir)/ws/test_web_services.py \
        --taxonomy-dir $(abs_top_srcdir)/data/ex-tax-1 \
        --synthesis-parent $(abs_top_srcdir)/data/ex-synth-par \
        --tests-parent $(abs_top_srcdir)/expectedws \
        --exe-dir $(abs_builddir)

