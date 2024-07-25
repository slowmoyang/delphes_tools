project(
  'delphes_tools',
  'cpp',
  license: 'MIT',
  version: '1.0.0',
)

#------------------------------------------------------------------------------
# dependencies
#------------------------------------------------------------------------------
# root
root_config = find_program('root-config', required: true)
root_compile_args = run_command(root_config, '--cflags', check: true).stdout().strip().split(' ')
root_link_args = run_command(root_config, '--ldflags', '--libs', check: true).stdout().strip().split(' ')
root_incdir = run_command(root_config, '--incdir', check: true).stdout().strip().split(' ')
root_dep = declare_dependency(
    compile_args: root_compile_args,
    link_args: root_link_args + ['-lEG'],
    include_directories: root_incdir
)

# delphes
delphes_config = find_program('delphes-config.py', required: true)
delphes_compile_args = run_command(delphes_config, '--cflags', check: true).stdout().strip().split(' ')
delphes_link_args = run_command(delphes_config, '--ldflags', '--libs', check: true).stdout().strip().split(' ')
delphes_incdir = run_command(delphes_config, '--incdir', check: true).stdout().strip().split(' ')
delphes_dep = declare_dependency(
  compile_args: delphes_compile_args,
  link_args: delphes_link_args,
  include_directories: delphes_incdir,
)

#------------------------------------------------------------------------------
# dependencies
#------------------------------------------------------------------------------
sources = [
  'src/genparticle.cc',
  'src/utils.cc',
  'src/weight.cc',
  'src/selection.cc',
]

delphes_tools_dep = declare_dependency(
  include_directories: include_directories('include'),
  sources: sources,
  dependencies: [
    root_dep,
    delphes_dep,
  ]
)

meson.override_dependency(meson.project_name(), delphes_tools_dep)
