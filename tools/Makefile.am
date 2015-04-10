LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
bin_PROGRAMS = otc-assert-invariants \
				otc-check-taxonomic-nodes \
				otc-detect-contested \
				otc-find-unsupported-nodes \
				otc-prune-taxonomy \
				otc-scaffolded-supertree \
				otc-count-leaves \
				otc-degree-distribution \
				otc-induced-subtree \
				otc-missing-splits \
				otc-polytomy-count \
				otc-taxon-conflict-report \
				otc-uncontested-decompose \
				otc-rf \
				otc-suppress-monotypic

otc_suppress_monotypic_SOURCES = suppressmonotypic.cpp
otc_suppress_monotypic_CPPFLAGS = $(AM_CPPFLAGS)

otc_assert_invariants_SOURCES = assertinvariants.cpp
otc_assert_invariants_CPPFLAGS = $(AM_CPPFLAGS)

otc_check_taxonomic_nodes_SOURCES = checktaxonomicnodes.cpp
otc_check_taxonomic_nodes_CPPFLAGS = $(AM_CPPFLAGS)

otc_detect_contested_SOURCES = detectcontested.cpp
otc_detect_contested_CPPFLAGS = $(AM_CPPFLAGS)

otc_uncontested_decompose_SOURCES = uncontesteddecompose.cpp
otc_uncontested_decompose_CPPFLAGS = $(AM_CPPFLAGS)

otc_find_unsupported_nodes_SOURCES = findunsupportednodes.cpp
otc_find_unsupported_nodes_CPPFLAGS = $(AM_CPPFLAGS)

otc_prune_taxonomy_SOURCES = prunetaxonomy.cpp
otc_prune_taxonomy_CPPFLAGS = $(AM_CPPFLAGS)

otc_scaffolded_supertree_SOURCES = scaffoldedsupertree.cpp
otc_scaffolded_supertree_CPPFLAGS = $(AM_CPPFLAGS)

otc_polytomy_count_SOURCES = polytomycount.cpp
otc_polytomy_count_CPPFLAGS = $(AM_CPPFLAGS)

otc_degree_distribution_SOURCES = degreedistribution.cpp
otc_degree_distribution_CPPFLAGS = $(AM_CPPFLAGS)

otc_count_leaves_SOURCES = countleaves.cpp
otc_count_leaves_CPPFLAGS = $(AM_CPPFLAGS)

otc_induced_subtree_SOURCES = inducedsubtree.cpp
otc_induced_subtree_CPPFLAGS = $(AM_CPPFLAGS)

otc_rf_SOURCES = rf.cpp
otc_rf_CPPFLAGS = $(AM_CPPFLAGS)

otc_missing_splits_SOURCES = missingsplits.cpp
otc_missing_splits_CPPFLAGS = $(AM_CPPFLAGS)

otc_taxon_conflict_report_SOURCES = taxonconflictreport.cpp
otc_taxon_conflict_report_CPPFLAGS = $(AM_CPPFLAGS)

check:
	python $(abs_top_srcdir)/tools/test_otc_tools.py $(abs_top_srcdir)/data $(abs_top_srcdir)/expected $(abs_builddir)
