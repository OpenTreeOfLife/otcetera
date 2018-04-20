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
  ['sampletips',  'sampletips'],
  ['displayedstats',  'displayed-stats'],
  ['subproblemstats',  'subproblem-stats'],
  ['nonterminaltoexemplars',   'nonterminals-to-exemplars'],
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
  ['solve-subproblem',  'solve-subproblem'],
  ['graft-solutions',  'graft-solutions'],
  ['unprune-solution',  'unprune-solution'],
  ['name-unnamed-nodes',  'name-unnamed-nodes'],
  ['annotate-synth',  'annotate-synth'],
  ['conflict-stats',  'conflict-stats'],
  ['taxonomy-parser',  'taxonomy-parser'],
  ['regrafttaxonomygenerator',  'regraft-taxonomy-generator'],
  ['relabel-tree',  'relabel-tree'],
  ['tree-tool',  'tree-tool'],
  ['broken-taxa',  'broken-taxa'],
  ['unprune-solution-and-name-unnamed-nodes', 'unprune-solution-and-name-unnamed-nodes']
  ]

foreach program : programs
  executable('otc-'+program[1], program[0] + '.cpp', dependencies: [boost, libotcetera, json], install_rpath: rpath, install: true)
endforeach

executable('otc-version-reporter', ['version-reporter.cpp',git_version_h], dependencies: [boost, libotcetera, json], install_rpath: rpath, install: true)
