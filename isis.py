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

def addparams( cmd, params ):
    if params:
        for name in params:
            cmd.append( '{}={}'.format(name, params[name])
    return cmd

def getkey(cube, grpname, keyword):
    '''Runs ISIS3 getkey'''
    # getkey grpname= ? keyword= ? from= ?
    print('Running getkey')
    path = os.environ['ISISROOT']+'/bin/getkey'
    result = subprocess.Popen([path, "from= "+ cube, "grpname= "+ grpname, "keyword= "+keyword], stdout=subprocess.PIPE).communicate()[0].strip()
    return result


def hi2isis(img, to=None, **keywords)
    '''Runs ISIS3 hi2isis'''
    print( 'Running hi2isis' )
    path = os.environ['ISISROOT']+'/bin/hi2isis'
    if to is None:
        to = os.path.splitext( img )[0] + '.hi2isis.cub'

    cmd = [path, 'from='+img, 'to='+to]
    r = subprocess.run( addparams(cmd, keywords), check=True )
    return

def histat( cub, to=None, **keywords )
    '''Runs ISIS3 histat'''
    print( 'Running histat' )
    path = os.environ['ISISROOT']+'/bin/histat'
    cmd = [path, 'from='+cub]
    cmd = addparams(cmd, keywords)
    if to is None:
        r = subprocess.run(cmd, capture_output=True, text=True check=True)
        return(r.stderr)
    else:
        cmd.append('to='+to)
        r = subprocess.run(cmd, check=True)
    return
