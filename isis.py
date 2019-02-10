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

import collections, csv, os, subprocess, sys

# These definitions and the use of env= in the subprocess.run calls allow us to
# run ISIS even though the shell that called this Python program may not be an
# ISIS-enabled shell.
isisroot = '/Users/rbeyer/.anaconda3/envs/isis3'
isis3data = '/Users/rbeyer/.anaconda3/envs/isis3/data'
isis_env = {'ISISROOT': isisroot, 
            'ISIS3DATA': isis3data, 
            'PATH': isisroot+'/bin/', 
            'HOME': os.environ['HOME']} # Otherwise ISIS tries to make a ./\$HOME dir


class Histogram:
    """A class to read and wrap the contents of the output of hist()"""

    def __init__(self, histfile):
        self.histfile = histfile
        (self.dictionary, 
         self.headers, 
         self.hist_list) = self.parsehist(histfile)
         # self.hist_list is a list of collections.namedtuple( 'HistRow', self.headers )
        
    def __str__(self):
        pass

    def __repr__(self):
        return (f'{self.__class__.__name__}(\'{self.histfile}\')')

    def __len__(self):
        return len(self.hist_list)

    def __getitem__(self,key):
        #if key is integer or slice object, look in self.hist_vals
        if isinstance(key, int):
            return self.hist_list[key]

        elif isinstance(key, slice):
            return self.hist_list[key]

        elif isinstance(key, str):
            return self.dictionary[key]

    def __iter__(self):
        return self.hist_list.__iter__()

    def __contains__(self, item):
        if item in self.dictionary: return True
        else: return False

    def keys(self):
        return self.dictionary.keys()

    def values(self):
        return self.dictionary.values()

    def parsehist(self, histfile):
        d = dict()
        headers = []
        hist_vals= []
        with open( histfile ) as f:
            for line in f:
                if ':' in line:
                    (k,v) = line.split(':')
                    if 'Cube' in k: d[k.strip()] = v.strip()
                    else:           d[k.strip()] = float(v.strip())
                    
                if line.startswith('DN,'): 
                    headers = line.strip().split(',')
                    HistRow = collections.namedtuple( 'HistRow', headers )
                    for row in map( HistRow._make, csv.reader(f, quoting=csv.QUOTE_NONNUMERIC) ):
                        hist_vals.append( row )

        return d, headers, hist_vals


#########################################################################
# Helper and wrapper functions for the ISIS commands.

def addparams( cmd, params ):
    '''Builds a list of strings from dictionary keys and values where the elements are "key=value".'''
    if params:
        for name in params:
            cmd.append( f'{name}={params[name]}' )
    return cmd

def _run_isis_program( cmd, capture=False ):
    '''Wrapper for subprocess.run()'''
    if capture: 
        return subprocess.run(cmd, env=isis_env, check=True, capture_output=True, text=True)
    else:
        return subprocess.run(cmd, env=isis_env, check=True)


def _build_isis_fn( fn_name ):
    '''This automatically builds a simple function to call an ISIS program, based on the name given.'''
    # Define the structure of the generic function, fn:
    def fn( fromcube, **keywords ):
        cmd = [fn_name, 'from='+fromcube]
        return( _run_isis_program( addparams(cmd, keywords) ) )
    fn.__name__ = fn_name
    fn.__doc__ = f'Runs ISIS3 {fn_name}'

    # Then add it, by name to the enclosing module.
    setattr( sys.modules[__name__], fn_name, fn)
    # Could have also used sys.modules['isis'] if I wanted to be specific.
       
# Now use the builder function to automatically create functions with these names:
_isis_programs = ['crop','hist','mask']
# Could also reach out and grab these names from the $ISISROOT/bin/ directory.
for p in _isis_programs:
    _build_isis_fn( p )


# Straight-up calls to ISIS programs below here, alphabetically arranged.
# These I want to call differently than the default builder provides.

def getkey(cube, grpname, keyword):
    '''Runs ISIS3 getkey'''
    cmd = ['getkey', 'from='+cube, 'grpname='+grpname, 'keyword='+keyword]
    return( _run_isis_program( cmd, capture=True ).stdout.strip() )
 
def hi2isis(img, to=None, **keywords):
    '''Runs ISIS3 hi2isis'''
    if to is None:
        to = os.path.splitext( img )[0] + '.hi2isis.cub'

    cmd = ['hi2isis', 'from='+img, 'to='+to]
    return( _run_isis_program( addparams(cmd, keywords) ) )

def histat(cub, to=None, **keywords):
    '''Runs ISIS3 histat'''
    cmd = ['histat', 'from='+cub]
    cmd = addparams(cmd, keywords)
    if to is None:
        return( _run_isis_program( cmd, capture=True ).stdout.strip() )
    else:
        cmd.append('to='+to)
        return( _run_isis_program( cmd ) )
