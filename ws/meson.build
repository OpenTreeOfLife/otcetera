

otc_tol_ws_sources = ['tolwsbooting.cpp',
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
     args: ['--taxonomy-dir',tax_dir,
     	    '--synthesis-parent',synth_dir,
	    	'--tests-parent',expectedws_dir,
	    	'--exe-dir',exe_dir,
			'--secs-to-recheck-pid-file=30']
    )

tax_dir = join_paths(meson.source_root(),'test/data/tnrs-taxonomy-induced')
synth_dir = join_paths(meson.source_root(),'test/data/tnrs-taxonomy-tree')
expectedws_dir = join_paths(meson.source_root(),'test/expectedws-induced')
test('tnrs web services test',
     test_web_services,
     timeout: 300,
     args: ['--taxonomy-dir', tax_dir,
     	    '--synthesis-parent', synth_dir,
	    	'--tests-parent', expectedws_dir,
	    	'--exe-dir', exe_dir,
			'--server-port=1986',
			'--secs-to-recheck-pid-file=30']
    )
