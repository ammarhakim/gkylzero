"""
Top-level paths
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('-p', type='string', help='Path to Gkyl dependency directory', dest='gkylDepsDir',
                   default=os.path.expandvars('$HOME/gkylsoft'))
    opt.add_option('--cxxflags', type='string', help='Compiler flags', dest='gkcxxflags',
                   default='-O3,-g,-Wall,-Wextra,-Wno-unused-variable,-Wno-unused-parameter,-std=c++17')
    opt.add_option('--cflags', type='string', help='Compiler flags', dest='gkcflags',
                   default='-O3,-g,-Wall,-Wextra,-Wno-unused-variable,-Wno-unused-parameter')
    opt.add_option('--debug', help='Debug flags', dest='gkdebug',
                   action='store_true', default=False)
    opt.add_option('--prefix', type='string', help='Install path', dest='prefix',
                   default=os.path.expandvars('$HOME/gkylsoft/gkylzero'))
    opt.add_option('--real-type', type='string', help='Real number type (float or double)', dest='real_type',
                   default="double")
    opt.add_option('--extra-link-libs', type='string', help='Extra libraries to link to', dest='extralibs',
                   default='')

@conf
def check_gkylzero(conf):
    conf.start_msg("Setting dependency path:")
    conf.end_msg(conf.options.gkylDepsDir)

    conf.start_msg("Setting prefix:")
    conf.end_msg(conf.options.prefix)
    conf.env.PREFIX = conf.options.prefix

    conf.env.append_value('CXXFLAGS', conf.options.gkcxxflags.split(','))
    conf.env.append_value('CFLAGS', conf.options.gkcflags.split(','))
    conf.env.append_value('LDFLAGS', conf.options.gkcflags.split(','))
    if conf.options.gkdebug:
      conf.env.append_value('CXXFLAGS', '-g')
      conf.env.append_value('CFLAGS', '-g')

    if conf.options.real_type == "float":
        conf.start_msg("Using real-type:")
        conf.end_msg("float")

        conf.env.append_value('CFLAGS', '-DGKYL_REAL_TYPE=1')
        conf.env.append_value('CXXFLAGS', '-DGKYL_REAL_TYPE=1')

    else:
        conf.start_msg("Using real-type:")
        conf.end_msg("double")

    conf.env.EXTRALIBS = ' '.join(conf.options.extralibs.split(','))
    
    return 1

def detect(conf):
    return detect_gkyl(conf)
