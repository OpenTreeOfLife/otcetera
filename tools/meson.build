git_version_sh = find_program('../scripts/git_version.sh')

git_version_h1 = vcs_tag(
  input : 'git_version.h.in',
  output: 'git_version.h.in.1',
  command: [ 'git', 'log', '-n1', '--date=format:%b %d %Y %k:%M:%S', '--pretty=format:[%D commit %h]  (%cd)' ],
  replace_string: '@git_message@',
  fallback: '',
)

git_version_h = vcs_tag(
  input : git_version_h1,
  output: 'git_version.h',
  command: [ 'git', 'rev-parse', 'HEAD'],
  replace_string: '@git_shaid@',
  fallback: '',
)

programs = [
  ['mungenames',  'munge-names'],
  ['sampletips',  'sample-tips'],
  ['displayedstats',  'displayed-stats'],
  ['subproblemstats',  'subproblem-stats'],
  ['nonterminaltoexemplars',   'nonterminals-to-exemplars'],
  ['moveextincthighertoavoidcontestingtaxa',   'move-extinct-higher-to-avoid-contesting-taxa'],
  ['setofids',  'set-of-ids'],
  ['suppressmonotypic',   'suppress-monotypic'],
  ['assertinvariants',   'assert-invariants'],
  ['detectcontested',  'detect-contested'],
  ['uncontesteddecompose',  'uncontested-decompose'],
  ['findresolution',  'find-resolution'],
  ['checksupertree',  'check-supertree'],
  ['prunetaxonomy',  'prune-taxonomy'],
  ['prunetosubtree',  'prune-to-subtree'],
  ['prunesynthtosubproblem',  'prune-synth-to-subproblem'],
  ['polytomycount',  'polytomy-count'],
  ['degreedistribution',  'degree-distribution'],
  ['countleaves',  'count-leaves'],
  ['inducedsubtree',  'induced-subtree'],
  ['distance',  'distance'],
  ['taxonconflictreport',  'taxon-conflict-report'],
  ['graft-solutions',  'graft-solutions'],
  ['unprune-solution',  'unprune-solution'],
  ['name-unnamed-nodes',  'name-unnamed-nodes'],
  ['annotate-synth',  'annotate-synth'],
  ['conflict-stats',  'conflict-stats'],
  ['source-taxonomy-patcher',  'source-taxonomy-patcher'],
  ['taxonomy-diff-maker',  'taxonomy-diff-maker'],
  ['taxonomy-parser',  'taxonomy-parser'],
  ['taxonomy-patcher',  'taxonomy-patcher'],
  ['regrafttaxonomygenerator',  'regraft-taxonomy-generator'],
  ['relabel-tree',  'relabel-tree'],
  ['tree-tool',  'tree-tool'],
  ['prune-clean',  'prune-clean'],
  ['broken-taxa',  'broken-taxa'],
  ['unprune-solution-and-name-unnamed-nodes', 'unprune-solution-and-name-unnamed-nodes'],
  ['checknamednodes', 'check-named-nodes'],
  ]

# we need restbed for this, indirectly.
if get_option('webservices')
  programs += [['tnrs-cli', 'tnrs-cli']]
endif

foreach program : programs
  executable('otc-'+program[1],
             program[0] + '.cpp',
             dependencies: [boost, libotcetera, json],
             install_rpath: rpath,
             install: true)
endforeach

solver_sources = ['solve-subproblem.cpp','solver/rsplit.cpp', 'solver/solution.cpp',
                  'solver/build.cpp', 'solver/rollback.cpp', 'solver/oracle.cpp',
                  'solver/combine.cpp', 'solver/names.cpp']

executable('otc-solve-subproblem',
           solver_sources,
           dependencies: [boost, libotcetera, json],
           install_rpath: rpath,
           install: true)

executable('otc-version-reporter', ['version-reporter.cpp',git_version_h], dependencies: [boost, libotcetera, json], install_rpath: rpath, install: true)



test_otc_tools = find_program(meson.project_source_root()/'tools/test_otc_tools.py')
data_dir =     meson.project_source_root()/'data'
expected_dir = meson.project_source_root()/'expected'
tools_dir =    meson.project_build_root()/'tools'

test('otc tools test',
     test_otc_tools,
     timeout: 300,
     workdir: tools_dir,
     args: [data_dir, expected_dir, tools_dir]
    )
