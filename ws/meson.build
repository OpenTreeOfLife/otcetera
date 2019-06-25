# We don't need ssl and curl unless librested itself was compiled with them.
#ssl  = dependency('openssl')
#curl = dependency('libcurl')
# Can we test if we need them to link with restbed?
threads = dependency('threads')

restbed_lib = cpp.find_library('restbed')
restbed_deps = [restbed_lib, threads] # ssl, curl?
restbed = declare_dependency(dependencies: restbed_deps)

# FIXME: we could put this back if we made enabling webservices a build-if-you-can tri-state option.
# if not restbed.found() or not ssl.found() or not curl.found()
#  restbed = disabler()
# endif

otc_tol_ws_sources = ['conflictws.cpp',
					  'taxonomyws.cpp',
	                  'tnrsws.cpp',
					  'tolws.cpp',
					  'tolwsadaptors.cpp',
					  'tolwsbooting.cpp',
					  'nexson/nexson.cpp',
					  'trees_to_serve.cpp',
					  ]
ws_inc = include_directories('.')
executable('otc-tol-ws', otc_tol_ws_sources,
	   dependencies: [boost, libotcetera, json, restbed],
	   include_directories: ws_inc,
	   install_rpath: rpath,
	   install: true)


test_web_services = find_program(join_paths(meson.source_root(),'ws/test_web_services.py'))

tax_dir = join_paths(meson.source_root(),'data/ex-tax-1')
synth_dir = join_paths(meson.source_root(),'data/ex-synth-par')
expectedws_dir = join_paths(meson.source_root(),'expectedws')
# this should be the dir of the otc-tol-ws binary before its installed.
exe_dir = join_paths(meson.build_root(),'ws')

test('web services test',
     test_web_services,
     timeout: 300,
#     workdir: exe_dir,
     args: ['--taxonomy-dir',tax_dir,
	    '--synthesis-parent',synth_dir,
	    '--tests-parent',expectedws_dir,
	    '--exe-dir',exe_dir]
    )
