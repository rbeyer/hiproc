#!/usr/bin/env python
"""This module contains calls to ISIS3 functions."""

# Copyright 2019, Ross A. Beyer (rbeyer@seti.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Someday there might be a real interface to ISIS3 through Python.  For now, we fake it:

import os, subprocess

# These definitions and the use of env= in the subprocess.run calls allow us to
# run ISIS even though the shell that called this Python program may not be an
# ISIS-enabled shell.
isisroot = '/Users/rbeyer/.anaconda3/envs/isis3'
isis3data = '/Users/rbeyer/.anaconda3/envs/isis3/data'
isis_env = {'ISISROOT': isisroot, 
            'ISIS3DATA': isis3data, 
            'PATH': isisroot+'/bin/', 
            'HOME': os.environ['HOME']} # Otherwise ISIS tries to make a ./\$HOME dir

def addparams( cmd, params ):
    if params:
        for name in params:
            cmd.append( f'{name}={params[name]}' )
    return cmd

def getkey(cube, grpname, keyword):
    '''Runs ISIS3 getkey'''
    cmd = ['getkey', 'from='+cube, 'grpname='+grpname, 'keyword='+keyword]
    r = subprocess.run( cmd, env=isis_env, check=True, capture_output=True, text=True )
    return( r.stdout.strip() )
 
def hi2isis(img, to=None, **keywords):
    '''Runs ISIS3 hi2isis'''
    if to is None:
        to = os.path.splitext( img )[0] + '.hi2isis.cub'

    cmd = ['hi2isis', 'from='+img, 'to='+to]
    subprocess.run( addparams(cmd, keywords), env=isis_env, check=True )
    return

def hist( cub, to=hist, **keywords):
    cmd = ['hist', 'from='+img, 'to='+to]
    subprocess.run( addparams(cmd, keywords), env=isis_env, check=True )
    return
 
def histat(cub, to=None, **keywords):
    '''Runs ISIS3 histat'''
    cmd = ['histat', 'from='+cub]
    cmd = addparams(cmd, keywords)
    if to is None:
        r = subprocess.run(cmd, env=isis_env, capture_output=True, text=True, check=True)
        return(r.stdout)
    else:
        cmd.append('to='+to)
        subprocess.run(cmd, env=isis_env, check=True)
    return
