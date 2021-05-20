"""
Detect MPI includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-mpi', help=('Enable parallel build'),
                   dest='enable_mpi', action='store_true',
                   default=True)
    opt.add_option('--disable-mpi', help=('Disable parallel build'),
                   dest='enable_mpi', action='store_false')
    opt.add_option('--mpi-inc-dir', type='string', help='Path to MPI includes', dest='mpiIncDir')
    opt.add_option('--mpi-lib-dir', type='string', help='Path to MPI libraries', dest='mpiLibDir')
    opt.add_option('--mpi-link-libs', type='string', help='List of MPI libraries to link',
                   dest='mpiLinkLibs', default='mpi')

@conf
def check_mpi(conf):
    opt = conf.options
    conf.env['MPI_FOUND'] = False
    
    if not conf.options.enable_mpi:
        return
    if conf.options.mpiIncDir:
        conf.env.INCLUDES_MPI = conf.options.mpiIncDir
    else:
        conf.env.INCLUDES_MPI = conf.options.gkylDepsDir+"/openmpi-3.1.2/include"

    if conf.options.mpiLibDir:
        conf.env.LIBPATH_MPI = conf.options.mpiLibDir
    else:
        conf.env.LIBPATH_MPI = conf.options.gkylDepsDir+"/openmpi-3.1.2/lib"

    libList = conf.options.mpiLinkLibs
    conf.env.LIB_MPI = libList.split(',')
    
    conf.start_msg('Checking for MPI')
    conf.check(header_name='mpi.h', features='cxx cxxprogram', use="MPI", mandatory=True)
    conf.end_msg("Found MPI")

    conf.env.append_value('CFLAGS', '-DHAVE_MPI')
    conf.env.append_value('CXXFLAGS', '-DHAVE_MPI')
    conf.env['MPI_FOUND'] = True
    return 1

def detect(conf):
    return detect_mpi(conf)
