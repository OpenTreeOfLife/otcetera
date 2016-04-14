LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
bin_PROGRAMS = testotcgreedyforest \
					testotcnewicktoken \
					testotctreefromnewick \
					testotctreeiter
LIBS = $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS)

testotcnewicktoken_SOURCES = test_otc_newicktoken.cpp
testotcnewicktoken_CPPFLAGS = $(AM_CPPFLAGS)

testotcgreedyforest_SOURCES = test_otc_greedyforest.cpp
testotcgreedyforest_CPPFLAGS = $(AM_CPPFLAGS)

testotctreefromnewick_SOURCES = test_otc_treefromnewick.cpp
testotctreefromnewick_CPPFLAGS = $(AM_CPPFLAGS)

testotctreeiter_SOURCES = test_otc_tree_iter.cpp
testotctreeiter_CPPFLAGS = $(AM_CPPFLAGS)

check:
	./testotcgreedyforest $(abs_top_srcdir)/data
	./testotcnewicktoken $(abs_top_srcdir)/data
	./testotctreefromnewick $(abs_top_srcdir)/data
	./testotctreeiter $(abs_top_srcdir)/data

