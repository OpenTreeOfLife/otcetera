SUBDIRS = \
	assert-invariants \
	check-taxonomic-nodes \
	detect-contested \
	find-unsupported-nodes \
	prune-singletons \
	simple-stats

check:
	python $(abs_top_srcdir)/tools/test_otc_tools.py $(abs_top_srcdir)/data $(abs_top_srcdir)/expected $(abs_builddir)

