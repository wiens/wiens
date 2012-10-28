VERSION='0.0.1'
APPNAME='demo'

top = '.'
out = 'built'


def options(opt):
        opt.load('compiler_cxx boost')

def configure(conf):
        conf.load('compiler_cxx boost')
        conf.check_boost('')
        conf.check_cxx(lib='lua', uselib_store='LUA')

def build(bld):
        mydir=bld.path.find_dir(top)
        bld.env.append_value('CXXFLAGS', ['-I' + mydir.abspath()+'/headers'])

        bld.program(source='phase_space_evolution_demo.cpp',
                    cxxflags='-Wall -pedantic -std=c++0x -O3 -g3',
                    target='phase_space_evolution_demo', use='BOOST LUA')
