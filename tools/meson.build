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
  'mungenames',
  'sampletips',
  'displayedstats',
  'subproblemstats',
  'nonterminaltoexemplars',
  'setofids',
  'suppressmonotypic',
  'assertinvariants',
  'detectcontested',
  'uncontesteddecompose',
  'findresolution',
  'checksupertree',
  'prunetaxonomy',
  'prunetosubtree',
  'prunesynthtosubproblem',
  'polytomycount',
  'degreedistribution',
  'countleaves',
  'inducedsubtree',
  'distance',
  'taxonconflictreport',
  'solve-subproblem',
  'graft-solutions',
  'unprune-solution',
  'name-unnamed-nodes',
  'annotate-synth',
  'conflict-stats',
  'taxonomy-parser',
  'regrafttaxonomygenerator',
  'relabel-tree',
  'tree-tool',
  'broken-taxa',
  'unprune-solution-and-name-unnamed-nodes'
  ]

foreach program : programs
  executable('otc-'+program, program + '.cpp', dependencies: [boost, libotcetera, json], include_directories: otc_inc, install: true)
endforeach

executable('otc-version-reporter', ['version-reporter.cpp',git_version_h], dependencies: [boost, libotcetera, json], include_directories: otc_inc, install: true)
