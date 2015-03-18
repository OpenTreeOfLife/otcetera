LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
noinst_PROGRAMS = testotcgreedyforest \
					testotcnewicktoken \
					testotctreefromnewick \
					testotctreeiter

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

