arpack_sources = [
  'ARPACK/_arpack_n_double_complex.c',
  'ARPACK/_arpack_n_double_complex.h',
  'ARPACK/_arpack_n_double.c',
  'ARPACK/_arpack_n_double.h',
  'ARPACK/_arpack_n_single_complex.c',
  'ARPACK/_arpack_n_single_complex.h',
  'ARPACK/_arpack_n_single.c',
  'ARPACK/_arpack_n_single.h',
  'ARPACK/_arpack_s_double.c',
  'ARPACK/_arpack_s_double.h',
  'ARPACK/_arpack_s_single.c',
  'ARPACK/_arpack_s_single.h',
  'ARPACK/_arpack.h'
]

arpack_lib = static_library('_arpack_lib',
  arpack_sources,
  fortran_args: [fortran_ignore_warnings, _suppress_all_warnings],
  include_directories: ['ARPACK/SRC'],
  override_options: ['b_lto=false'],
  gnu_symbol_visibility: 'hidden',
)

arpack_module = custom_target('arpack_module',
  output: ['_arpackmodule.c', '_arpack-f2pywrappers.f'],
  input: 'arpack.pyf.src',
  command: [generate_f2pymod, '@INPUT@', '-o', '@OUTDIR@'] + f2py_freethreading_arg,
)

_arpack = py3.extension_module('_arpack',
  arpack_module,
  link_with: arpack_lib,
  include_directories: ['../../../../_lib/src'],
  dependencies: [np_dep, lapack_dep, blas_dep],
  link_args: version_link_args,
  install: true
)


py3.install_sources([
    '__init__.py',
    'arpack.py',
  ],
  subdir: 'scipy/sparse/linalg/_eigen/arpack'
)

subdir('tests')
