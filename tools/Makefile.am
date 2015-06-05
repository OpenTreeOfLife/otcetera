LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc
bin_PROGRAMS = otc-assert-invariants \
				otc-count-leaves \
				otc-degree-distribution \
				otc-detect-contested \
				otc-displayed-stats \
				otc-distance \
				otc-check-supertree \
				otc-find-resolution \
				otc-induced-subtree \
				otc-nonterminals-to-exemplars \
				otc-polytomy-count \
				otc-prune-synth-to-subproblem \
				otc-prune-taxonomy \
				otc-prune-to-subtree \
				otc-scaffolded-supertree \
				otc-set-of-ids \
				otc-subproblem-stats \
				otc-suppress-monotypic \
				otc-taxon-conflict-report \
				otc-uncontested-decompose

otc_displayed_stats_SOURCES = displayedstats.cpp
otc_displayed_stats_CPPFLAGS = $(AM_CPPFLAGS)


otc_subproblem_stats_SOURCES = subproblemstats.cpp
otc_subproblem_stats_CPPFLAGS = $(AM_CPPFLAGS)

otc_nonterminals_to_exemplars_SOURCES = nonterminaltoexemplars.cpp
otc_nonterminals_to_exemplars_CPPFLAGS = $(AM_CPPFLAGS)

otc_set_of_ids_SOURCES = setofids.cpp
otc_set_of_ids_CPPFLAGS = $(AM_CPPFLAGS)

otc_suppress_monotypic_SOURCES = suppressmonotypic.cpp
otc_suppress_monotypic_CPPFLAGS = $(AM_CPPFLAGS)

otc_assert_invariants_SOURCES = assertinvariants.cpp
otc_assert_invariants_CPPFLAGS = $(AM_CPPFLAGS)

otc_detect_contested_SOURCES = detectcontested.cpp
otc_detect_contested_CPPFLAGS = $(AM_CPPFLAGS)

otc_uncontested_decompose_SOURCES = uncontesteddecompose.cpp
otc_uncontested_decompose_CPPFLAGS = $(AM_CPPFLAGS)

otc_find_resolution_SOURCES = findresolution.cpp
otc_find_resolution_CPPFLAGS = $(AM_CPPFLAGS)

otc_check_supertree_SOURCES = checksupertree.cpp
otc_check_supertree_CPPFLAGS = $(AM_CPPFLAGS)

otc_prune_taxonomy_SOURCES = prunetaxonomy.cpp
otc_prune_taxonomy_CPPFLAGS = $(AM_CPPFLAGS)

otc_prune_to_subtree_SOURCES = prunetosubtree.cpp
otc_prune_to_subtree_CPPFLAGS = $(AM_CPPFLAGS)

otc_prune_synth_to_subproblem_SOURCES = prunesynthtosubproblem.cpp
otc_prune_synth_to_subproblem_CPPFLAGS = $(AM_CPPFLAGS)

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

otc_distance_SOURCES = distance.cpp
otc_distance_CPPFLAGS = $(AM_CPPFLAGS)

otc_taxon_conflict_report_SOURCES = taxonconflictreport.cpp
otc_taxon_conflict_report_CPPFLAGS = $(AM_CPPFLAGS)

check:
	python $(abs_top_srcdir)/tools/test_otc_tools.py $(abs_top_srcdir)/data $(abs_top_srcdir)/expected $(abs_builddir)
