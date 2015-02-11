LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
noinst_PROGRAMS = testotcnewicktoken
testotcnewicktoken_SOURCES = test_otc_newicktoken.cpp
testotcnewicktoken_CPPFLAGS = $(AM_CPPFLAGS)

check:
	./testotcnewicktoken $(abs_top_srcdir)/data

