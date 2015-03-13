LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
bin_PROGRAMS = otcassertinvariants \
				otcchecktaxonomicnodes \
				otcdetectcontested \
				otcfindunsupportednodes \
				otcprunetaxonomy \
				otcremaptodeepestunlisted \
				otcscaffoldedsupertree \
				otccountleaves \
				otcdegreedistribution \
				otcinducedsubtree \
				otcmissingsplits \
				otcpolytomycount \
				otctaxonconflictreport \
				otcrf

otcassertinvariants_SOURCES = assertinvariants.cpp
otcassertinvariants_CPPFLAGS = $(AM_CPPFLAGS)

otcchecktaxonomicnodes_SOURCES = checktaxonomicnodes.cpp
otcchecktaxonomicnodes_CPPFLAGS = $(AM_CPPFLAGS)

otcdetectcontested_SOURCES = detectcontested.cpp
otcdetectcontested_CPPFLAGS = $(AM_CPPFLAGS)

otcfindunsupportednodes_SOURCES = findunsupportednodes.cpp
otcfindunsupportednodes_CPPFLAGS = $(AM_CPPFLAGS)

otcprunetaxonomy_SOURCES = prunetaxonomy.cpp
otcprunetaxonomy_CPPFLAGS = $(AM_CPPFLAGS)

otcremaptodeepestunlisted_SOURCES = remaptodeepestunlisted.cpp
otcremaptodeepestunlisted_CPPFLAGS = $(AM_CPPFLAGS)

otcscaffoldedsupertree_SOURCES = scaffoldedsupertree.cpp
otcscaffoldedsupertree_CPPFLAGS = $(AM_CPPFLAGS)

otcpolytomycount_SOURCES = polytomycount.cpp
otcpolytomycount_CPPFLAGS = $(AM_CPPFLAGS)

otcdegreedistribution_SOURCES = degreedistribution.cpp
otcdegreedistribution_CPPFLAGS = $(AM_CPPFLAGS)

otccountleaves_SOURCES = countleaves.cpp
otccountleaves_CPPFLAGS = $(AM_CPPFLAGS)

otcinducedsubtree_SOURCES = inducedsubtree.cpp
otcinducedsubtree_CPPFLAGS = $(AM_CPPFLAGS)

otcrf_SOURCES = rf.cpp
otcrf_CPPFLAGS = $(AM_CPPFLAGS)

otcmissingsplits_SOURCES = missingsplits.cpp
otcmissingsplits_CPPFLAGS = $(AM_CPPFLAGS)

otctaxonconflictreport_SOURCES = taxonconflictreport.cpp
otctaxonconflictreport_CPPFLAGS = $(AM_CPPFLAGS)

check:
	python $(abs_top_srcdir)/tools/test_otc_tools.py $(abs_top_srcdir)/data $(abs_top_srcdir)/expected $(abs_builddir)
