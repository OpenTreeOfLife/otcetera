

otc_tol_ws_sources = ['tolwsbooting.cpp']

ws_inc = include_directories('.')
executable('otc-tol-ws', otc_tol_ws_sources,
	   dependencies: [boost, libotcetera, json, restbed],
	   include_directories: ws_inc,
	   install_rpath: rpath,
	   install: true)


test_web_services = find_program(join_paths(meson.source_root(),'ws/test_web_services.py'))

# this should be the dir of the otc-tol-ws binary before its installed.
exe_dir = join_paths(meson.build_root(),'ws')

synth_setups = join_paths(meson.source_root(),'data/synth-setups')

synth_setup_1 = join_paths(synth_setups,'synth-1')

tax_dir1        = join_paths(synth_setup_1, 'taxonomy')
synth_dir1      = join_paths(synth_setup_1, 'synth-par')
expectedws_dir1 = join_paths(synth_setup_1, 'expectedws')

test('web services test (1)',
     test_web_services,
     timeout: 300,
     args: ['--taxonomy-dir',tax_dir1,
            '--synthesis-parent',synth_dir1,
            '--tests-parent',expectedws_dir1,
            '--exe-dir',exe_dir,
            '--secs-to-recheck-pid-file=30']
    )

synth_setup_2 = join_paths(synth_setups,'synth-2')

tax_dir2        = join_paths(synth_setup_2, 'taxonomy')
synth_dir2      = join_paths(synth_setup_2, 'synth-par')
expectedws_dir2 = join_paths(synth_setup_2, 'expectedws')

test('web services test (2)',
     test_web_services,
     timeout: 300,
     args: ['--taxonomy-dir',tax_dir2,
            '--synthesis-parent',synth_dir2,
            '--tests-parent',expectedws_dir2,
            '--exe-dir',exe_dir,
            '--secs-to-recheck-pid-file=30']
    )

tax_dir = join_paths(meson.source_root(),'test/tnrs/data/taxonomy-induced')
synth_dir = join_paths(meson.source_root(),'test/tnrs/data/taxonomy-tree')
expectedws_dir = join_paths(meson.source_root(),'test/tnrs/expectedws-induced')
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
