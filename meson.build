project(
  'josephson',
  'c',
  meson_version : '>= 1.3.0',
  version : '0.1',
  default_options : ['warning_level=3'],
)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

libj_dep = library('libj', 'libj.c', dependencies: [m_dep])

dependencies = [
  m_dep,
]

links = [
  libj_dep,
]

exe = executable(
  'josephson',
  'josephson.c',
  dependencies : dependencies,
  link_with: links,
  install : true,
)

test('basic', exe)
