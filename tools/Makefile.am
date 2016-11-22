GIT_VERSION_CMD = $(SHELL) $(top_srcdir)/scripts/git_version.sh
git_version.stamp:
	@if test -f "$(srcdir)/git_version.h"; then \
		if test -f "git_version.h"; then :; \
		else \
			cp "$(srcdir)/git_version.h" "git_version.h"; \
		fi; \
	fi
	$(GIT_VERSION_CMD) -k -s $(top_srcdir) -o git_version.h
	@if test -s "$(srcdir)/git_version.h"; then \
		if cmp "$(srcdir)/git_version.h" "git_version.h"; then :; \
		else \
			echo "Error: $(srcdir)/git_version.h and git_version.h differ."; \
			echo "       You probably want to remove the former."; \
			exit 1; \
		fi; \
	fi

dist-hook: git_version.stamp
	if test -f "git_version.h"; then \
		$(SED) -e 's|^#undef GIT_IS_DIST.*|#define GIT_IS_DIST 1|' \
			"git_version.h" > "$(distdir)/git_version.h"; \
	fi

LDADD       = @top_builddir@/otc/libotcetera.la
AM_CPPFLAGS = -I@top_srcdir@/otc $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS)
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
				otc-set-of-ids \
				otc-subproblem-stats \
				otc-suppress-monotypic \
				otc-taxon-conflict-report \
				otc-uncontested-decompose \
				otc-solve-subproblem \
				otc-graft-solutions \
				otc-unprune-solution \
				otc-name-unnamed-nodes \
				otc-unprune-solution-and-name-unnamed-nodes \
				otc-annotate-synth \
				otc-taxonomy-parser \
				otc-relabel-tree \
				otc-tree-tool \
				otc-version-reporter \
				otc-munge-names \
				otc-conflict-stats



LIBS = $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS) 

otc_munge_names_SOURCES = mungenames.cpp
otc_munge_names_CPPFLAGS = $(AM_CPPFLAGS)


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

otc_solve_subproblem_SOURCES = solve-subproblem.cpp
otc_solve_subproblem_CPPFLAGS = $(AM_CPPFLAGS)

otc_graft_solutions_SOURCES = graft-solutions.cpp
otc_graft_solutions_CPPFLAGS = $(AM_CPPFLAGS)

otc_unprune_solution_SOURCES = unprune-solution.cpp
otc_unprune_solution_CPPFLAGS = $(AM_CPPFLAGS)

otc_name_unnamed_nodes_SOURCES = name-unnamed-nodes.cpp
otc_name_unnamed_nodes_CPPFLAGS = $(AM_CPPFLAGS)

otc_annotate_synth_SOURCES = annotate-synth.cpp
otc_annotate_synth_CPPFLAGS = $(AM_CPPFLAGS)

otc_conflict_stats_SOURCES = conflict-stats.cpp
otc_conflict_stats_CPPFLAGS = $(AM_CPPFLAGS)

otc_taxonomy_parser_SOURCES = taxonomy-parser.cpp
otc_taxonomy_parser_CPPFLAGS = $(AM_CPPFLAGS)

otc_relabel_tree_SOURCES = relabel-tree.cpp
otc_relabel_tree_CPPFLAGS = $(AM_CPPFLAGS)

otc_tree_tool_SOURCES = tree-tool.cpp
otc_tree_tool_CPPFLAGS = $(AM_CPPFLAGS)

otc_unprune_solution_and_name_unnamed_nodes_SOURCES = unprune-solution-and-name-unnamed-nodes.cpp
otc_unprune_solution_and_name_unnamed_nodes_CPPFLAGS = $(AM_CPPFLAGS)

otc_version_reporter_SOURCES = version-reporter.cpp
otc_version_reporter_CPPFLAGS = $(AM_CPPFLAGS)

nodist_otc_version_reporter_SOURCES = git_version.h

# always "rebuild" these
BUILT_SOURCES = version-reporter.cpp git_version.stamp

check:
	python $(abs_top_srcdir)/tools/test_otc_tools.py $(abs_top_srcdir)/data $(abs_top_srcdir)/expected $(abs_builddir)

