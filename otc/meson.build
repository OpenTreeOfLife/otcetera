libotcetera_sources = [
  'config_file.cpp',
  'ctrie/context_ctrie_db.cpp',
  'ctrie/str_utils.cpp',
  'ctrie/ctrie_db.cpp',
  'ctrie/ctrie_node.cpp',
  'ctrie/ctrie.cpp',
  'ctrie/search_impl.cpp',
  'easylogging++.cc',
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

libotcetera_ws_sources = [
  'ws/conflictws.cpp',
  'ws/taxonomyws.cpp',
  'ws/tnrsws.cpp',
  'ws/tolws.cpp',
  'ws/find_node.cpp',
  'ws/trees_to_serve.cpp', 
  'ws/nexson/nexson.cpp',
]

if get_option('webservices')
  libotcetera_sources += libotcetera_ws_sources
endif

otc_inc = include_directories('..')
otc_lib = shared_library('otcetera',
                         libotcetera_sources,
                         include_directories: otc_inc,
                         dependencies: [boost, restbed],
                         install: true)

libotcetera = declare_dependency(include_directories: otc_inc, link_with: otc_lib, dependencies: [restbed])

