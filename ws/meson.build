ssl  = dependency('openssl', required: false)
curl = dependency('libcurl', required: false)
threads = dependency('threads', required: false)

restbed = cpp.find_library('restbed', dependencies: [ssl, curl, threads], required: false)
if not restbed.found() or not ssl.found() or not curl.found()
  restbed = disabler()
endif

otc_tol_ws_sources = ['tolws.cpp', 'tolwsadaptors.cpp', 'tolwsbooting.cpp', 'nexson/nexson.cpp']
ws_inc = include_directories('.')
executable('otc-tol-ws', otc_tol_ws_sources,
	   dependencies: [boost, libotcetera, restbed, json, ssl, curl, threads],
	   include_directories: ws_inc,
	   install_rpath: rpath,
	   install: true)

#check:
#	python $(abs_top_srcdir)/ws/test_web_services.py \
#        --taxonomy-dir $(abs_top_srcdir)/data/ex-tax-1 \
#        --synthesis-parent $(abs_top_srcdir)/data/ex-synth-par \
#        --tests-parent $(abs_top_srcdir)/expectedws \
#        --exe-dir $(abs_builddir)

