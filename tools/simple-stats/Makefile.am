LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
bin_PROGRAMS = otccountleaves \
				otcdegreedistribution \
				otcinducedsubtree \
				otcpolytomycount 
otcpolytomycount_SOURCES = polytomycount.cpp
otcpolytomycount_CPPFLAGS = $(AM_CPPFLAGS)

otcdegreedistribution_SOURCES = degreedistribution.cpp
otcdegreedistribution_CPPFLAGS = $(AM_CPPFLAGS)

otccountleaves_SOURCES = countleaves.cpp
otccountleaves_CPPFLAGS = $(AM_CPPFLAGS)

otcinducedsubtree_SOURCES = inducedsubtree.cpp
otcinducedsubtree_CPPFLAGS = $(AM_CPPFLAGS)
