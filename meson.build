project('otcetera', ['cpp'],
	version: '0.0.01-dev',
	default_options: [
	  'buildtype=release',
          'b_ndebug=if-release',
	  'cpp_std=c++17',
          'warning_level=3'
	],
	meson_version: '>=0.53')

cpp = meson.get_compiler('cpp')

project_cpp_args = ['-Wno-sign-compare','-Woverloaded-virtual','-Wstrict-aliasing','-Wno-unknown-pragmas','-Wno-maybe-uninitialized','-Wno-expansion-to-defined']

add_project_arguments(cpp.get_supported_arguments(project_cpp_args), language: 'cpp')

# 2. Write a 'config.h'
conf_data = configuration_data()
conf_data.set_quoted('PACKAGE_VERSION', meson.project_version())
conf_data.set_quoted('_ARCH_', host_machine.system()+' ' + host_machine.cpu_family())
conf_data.set_quoted('_COMPILER_', cpp.get_id() + ' ' + cpp.version()+' ' + host_machine.cpu_family())
configure_file(output : 'config.h', configuration : conf_data)

# Do we need this?
# AC_PROG_LN_S

boost = dependency('boost', modules : ['program_options','system','filesystem'], version: '>=1.54')

json = declare_dependency(include_directories: include_directories('otc'))

if get_option('webservices')
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
else
  restbed = dependency('', required: false)
endif

rpath = '$ORIGIN/../'+get_option('libdir')

# build libotcetera
subdir('otc')
# build tools
subdir('tools')
# build web services
if get_option('webservices')
  subdir('ws')
endif

# check if python subprocess module is installed
if true
  subdir('test')  
endif

summary({'host': host_machine.system()
        }, section: 'Architecture')

summary({'webservices': get_option('webservices'),
        },section: 'Configuration')

summary({'prefix': get_option('prefix'),
        },section: 'Directories')

assertions_off = get_option('b_ndebug') == 'true' or (get_option('b_ndebug') == 'if-release' and
                                                      (get_option('buildtype') == 'release' or
                                                       get_option('buildtype') == 'plain')
                                                     )
assertions_enabled = 'enabled'
if assertions_off
  assertions_enabled = 'disabled'
endif

summary({'optimization': get_option('optimization'),
         'debug': get_option('debug'),
         'assertions': assertions_enabled,
        },section: 'Build')
