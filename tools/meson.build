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
  ['annotate-synth', 'annotate-synth'],
  ['assertinvariants', 'assert-invariants'],
  ['broken-taxa', 'broken-taxa'],
  ['checksupertree', 'check-supertree'],
  ['conflict-stats', 'conflict-stats'],
  ['contrast-astral-runs', 'contrast-astral-runs'],
  ['countleaves', 'count-leaves'],
  ['degreedistribution', 'degree-distribution'],
  ['detectcontested', 'detect-contested'],
  ['displayedstats', 'displayed-stats'],
  ['distance', 'distance'],
  ['explain-phylo-diffs', 'explain-phylo-diffs'],
  ['findresolution', 'find-resolution'],
  ['graft-solutions', 'graft-solutions'],
  ['inducedsubtree', 'induced-subtree'],
  ['moveextincthighertoavoidcontestingtaxa', 'move-extinct-higher-to-avoid-contesting-taxa'],
  ['mungenames', 'munge-names'],
  ['name-unnamed-nodes', 'name-unnamed-nodes'],
  ['nonterminaltoexemplars', 'nonterminals-to-exemplars'],
  ['polytomycount', 'polytomy-count'],
  ['prune-clean', 'prune-clean'],
  ['prunesynthtosubproblem', 'prune-synth-to-subproblem'],
  ['prunetaxonomy', 'prune-taxonomy'],
  ['prunetosubtree', 'prune-to-subtree'],
  ['quartet-distances', 'quartet-distances'],
  ['regrafttaxonomygenerator', 'regraft-taxonomy-generator'],
  ['relabel-tree', 'relabel-tree'],
  ['sampletips', 'sample-tips'],
  ['setofids', 'set-of-ids'],
  ['solve-subproblem', 'solve-subproblem'],
  ['subproblemstats', 'subproblem-stats'],
  ['suppressmonotypic', 'suppress-monotypic'],
  ['taxonconflictreport', 'taxon-conflict-report'],
  ['taxonomy-parser', 'taxonomy-parser'],
  ['tree-tool', 'tree-tool'],
  ['triple-distances', 'triple-distances'],
  ['uncontesteddecompose', 'uncontested-decompose'],
  ['unprune-solution', 'unprune-solution'],
  ['unprune-solution-and-name-unnamed-nodes', 'unprune-solution-and-name-unnamed-nodes'],
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

executable('otc-version-reporter', ['version-reporter.cpp',git_version_h], dependencies: [boost, libotcetera, json], install_rpath: rpath, install: true)



test_otc_tools = find_program(meson.source_root()/'tools/test_otc_tools.py')
data_dir =     meson.source_root()/'data'
expected_dir = meson.source_root()/'expected'
tools_dir =    meson.build_root()/'tools'

test('otc tools test',
     test_otc_tools,
     timeout: 300,
     workdir: tools_dir,
     args: [data_dir, expected_dir, tools_dir]
    )
