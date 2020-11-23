# -*- python -*-
# Top-level build-script for GkylZero
##    _______     ___
## + 6 @ |||| # P ||| +

import datetime
import os
import platform
import sys

APPNAME = 'gkylzero'
VER = "0.1"

now = datetime.datetime.now()
VERSION = VER + "-"+now.strftime("%Y-%m-%d")

top = '.'
out = 'build'

# extra flags to pass to linker
EXTRA_LINK_FLAGS = []

from waflib import TaskGen
from waflib.Options import options

def options(opt):
    opt.load('compiler_c compiler_cxx') 
    opt.load('gkylzero',
             tooldir='waf_tools')

def configure(conf):
    r"""Configure Gkyl build"""

    # load tools
    conf.load('compiler_c compiler_cxx')
    conf.check_gkylzero()

    # standard install location for dependencies
    gkydepsDir = os.path.expandvars('$HOME/gkylsoft')

    # add current build directory to pick up config header
    conf.env.append_value('INCLUDES', ['.'])
    
    # load options for math and dynamic library
    conf.env.LIB_M = ['m']
    conf.env.LIB_DL = ['dl']

    # write out configuration info into header
    conf.write_config_header('gkylconfig.h')


from waflib import Task
class GitTip(Task.Task):
    always_run = True # need to force running every time
    run_str = r'echo \#define GKYL_GIT_CHANGESET  \"`git describe --abbrev=12 --always --dirty=+`\" > ${TGT}'

def build(bld):
    gitTip = GitTip(env=bld.env)
    gitTip.set_outputs(bld.path.find_or_declare('gkylgittip.h'))
    bld.add_to_group(gitTip)

    if bld.jobs > 16:
      bld.jobs = 16
    
    # recurse down directories and build C code
    #bld.recurse("Comm")

    # build executable
    buildExec(bld)

def buildExec(bld):
    pass

def dist(ctx):
    ctx.algo = "zip" # use ZIP instead of tar.bz2
    ctx.excl = " **/.waf* **/*~ **/*.pyc **/*.swp **/.lock-w* build"
