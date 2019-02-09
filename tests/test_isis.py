#!/usr/bin/env python
"""This module has tests for the ISIS3 utility functions."""

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

import unittest, os, sys
sys.path.append('../')
from isis import *

class TestParams(unittest.TestCase):
    def test_addparams(self):
        t = ( 'isisprogram','from=foo.cub' )
        p = { 'to': 'to.cub', 'check': False, 'value': 3.0 }
        c = list(t)
        truth = list(t)
        truth.extend( ['to=to.cub', 'check=False', 'value=3.0'] )

        c = addparams( c, p )
        self.assertEqual( truth, c )

@unittest.skip('Takes a while to run hi2isis.')
class Test_hi2isis(unittest.TestCase):
    def setUp(self):
        self.img = 'test_hi2isis.img'

    def test_hi2isis_with_to(self):
        tocube = 'test_hi2isis.cub'
        hi2isis( self.img, tocube )
        self.assertTrue( os.path.isfile(tocube) )
        os.remove( tocube )

    def test_hi2isis_without_to(self):
        tocube = 'test_hi2isis.hi2isis.cub'
        hi2isis( self.img )
        self.assertTrue( os.path.isfile(tocube) )
        os.remove( tocube )

    def tearDown(self):
        os.remove( 'print.prt' )

@unittest.skip('Takes a while to run hi2isis.')
class Test_getkey(unittest.TestCase):
    def setUp(self):
        self.img = 'test_hi2isis.img'
        self.cub = 'test_getkey.cub'
        hi2isis( self.img, self.cub )
    
    def tearDown(self):
        os.remove( self.cub )
        os.remove( 'print.prt' )

    def test_getkey(self):
        truth = 'HIRISE'
        key = getkey( self.cub, 'Instrument', 'InstrumentId' )
        self.assertEqual( truth, key )

    def test_getkey_fail(self):
        # Pixels doesn't have InstrumentId, should fail
        self.assertRaises( subprocess.CalledProcessError, getkey, self.cub, 'Pixels', 'InstrumentId' )

@unittest.skip('Takes a while to run hi2isis.')
class Test_histat(unittest.TestCase):
    def setUp(self):
        self.img = 'test_hi2isis.img'
        self.cub = 'test_histat.cub'
        hi2isis( self.img, self.cub )
    
    def tearDown(self):
        os.remove( self.cub )
        os.remove( 'print.prt' )

    def test_histat_with_to(self):
        tofile = self.cub+'.histat'
        histat( self.cub, tofile )
        self.assertTrue( os.path.isfile(tofile) )
        os.remove( tofile )

    def test_histat_without_to(self):
        s = histat( self.cub )
        self.assertTrue( s.startswith('Group = IMAGE_POSTRAMP') )
