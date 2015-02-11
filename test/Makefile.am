LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
noinst_PROGRAMS = testotcnewicktoken
testotcnewicktoken_SOURCES = test_otc_newicktoken.cpp
testotcnewicktoken_CPPFLAGS = $(AM_CPPFLAGS)

#installcheck-local:
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py -x $(bindir)/NEXUSnormalizer $(srcdir)/funkyValidIn $(srcdir)/funkyValidOut
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py -e $(bindir)/NEXUSnormalizer $(srcdir)/ExternalValidIn $(srcdir)/ExternalValidOut
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py -i -e $(bindir)/NEXUSnormalizer $(srcdir)/ExternalInvalidIn
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py -a $(bindir)/NEXUSnormalizer
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/NEXUSnormalizer $(top_srcdir)/data/characters.nex $(srcdir)/data
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/NEXUSnormalizer $(top_srcdir)/data/sample.tre $(srcdir)/data
#	$(PYTHON) $(srcdir)/roundTripNCLTest.py $(bindir)/NEXUSnormalizer $(srcdir)/NTSValidIn $(srcdir)/NTSValidOut
