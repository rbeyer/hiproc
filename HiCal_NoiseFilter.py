#!/usr/bin/env python
"""Perform NoiseFilter processing."""

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


    # reads things from a conf file, pvl
    # sets up a bunch of filenames
    # isis.hist
    # GetHistVal() - reads info from histogram file
    # isis.cubenorm
    # read the resulting cubenorm.tab file
    # does some math to filter that info, then writes to a textfile
    # isis.cubenorm
    # isis.lowpass, highpass, algebra, noisefilter a bunch of times


if __name__ == "__main__":
    main()
