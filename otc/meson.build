libotcetera_sources = [
  'config_file.cpp',
  'ctrie/context_ctrie_db.cpp',
  'ctrie/str_utils.cpp',
  'embedded_tree.cpp',
  'forest.cpp',
  'ftree.cpp',
  'greedy_forest.cpp',
  'newick.cpp',
  'node_embedding.cpp',
  'otcetera.cpp',
  'otcli.cpp',
  'supertree_util.cpp',
  'taxonomy/taxonomy.cpp',
  'taxonomy/flags.cpp',
  'test_harness.cpp',
  'tnrs/nomenclature.cpp',
  'tnrs/context.cpp',
  'tree.cpp',
  'util.cpp',
  'write_dot.cpp',
  ]

otc_inc = include_directories('..')
otc_lib = shared_library('otcetera', libotcetera_sources, include_directories: otc_inc, dependencies: boost, install: true)

libotcetera = declare_dependency(include_directories: otc_inc, link_with: otc_lib)

