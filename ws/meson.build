

otc_tol_ws_sources = ['tolwsbooting.cpp']

ws_inc = include_directories('.')
executable('otc-tol-ws', otc_tol_ws_sources,
	   dependencies: [boost, libotcetera, json, restbed],
	   include_directories: ws_inc,
	   install_rpath: rpath,
	   install: true)


test_web_services = find_program( meson.project_source_root()/'ws/test_web_services.py' )

# this should be the dir of the otc-tol-ws binary before its installed.
exe_dir = meson.project_build_root()/'ws'

synth_setups = meson.project_source_root()/'data/synth-setups'

synth_setup_1 = synth_setups/'synth-1'

tax_dir1        = synth_setup_1/'taxonomy'
synth_dir1      = synth_setup_1/'synth-par'
expectedws_dir1 = synth_setup_1/'expectedws'

test('web services test (1)',
     test_web_services,
     timeout: 300,
     args: ['--taxonomy-dir',tax_dir1,
            '--synthesis-parent',synth_dir1,
            '--tests-parent',expectedws_dir1,
            '--exe-dir',exe_dir,
	    '--server-port=1987',
            '--secs-to-recheck-pid-file=30']
    )

synth_setup_2 = synth_setups/'synth-2'

tax_dir2        = synth_setup_2/'taxonomy'
synth_dir2      = synth_setup_2/'synth-par'
expectedws_dir2 = synth_setup_2/'expectedws'

test('web services test (2)',
     test_web_services,
     timeout: 300,
     args: ['--taxonomy-dir',tax_dir2,
            '--synthesis-parent',synth_dir2,
            '--tests-parent',expectedws_dir2,
            '--exe-dir',exe_dir,
	    '--server-port=1988',
            '--secs-to-recheck-pid-file=30']
    )

tax_dir = meson.project_source_root()/'test/tnrs/data/taxonomy-induced'
synth_dir = meson.project_source_root()/'test/tnrs/data/taxonomy-tree'
expectedws_dir = meson.project_source_root()/'test/tnrs/expectedws-induced'
test('tnrs web services test',
     test_web_services,
     timeout: 300,
     args: ['--taxonomy-dir', tax_dir,
     	    '--synthesis-parent', synth_dir,
	    '--tests-parent', expectedws_dir,
	    '--exe-dir', exe_dir,
	    '--server-port=1989',
	    '--secs-to-recheck-pid-file=30']
    )
