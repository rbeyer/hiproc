.. highlight:: shell

============
Installation
============


.. Stable release
   --------------

To install hiproc, run this command in your terminal:

    .. code-block:: console

        $ pip install hiproc

This is the preferred method to install hiproc, as it will always
install the most recent stable release.

However, hiproc depends on the Python gdal implementation, which has a C library
requirement, and may not install properly if it is not present, and doesn't always
smoothly install via pip alone.  If you have issues with hiproc installation, you
may want to try installing the Python gdal separately or via conda, and then attempt
to pip install hiproc again.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

Eventually I'll get a conda install going.

This package depends on `kalasiris`_ which is just a wrapper around
`ISIS`_, which has some complicated interactions with other
environments.  You will need to install ISIS separately in another
conda environment, and then review the strategies detailed in the
kalasiris page on `ISIS Interaction`_ to ensure that you can get
everything working.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/
.. _kalasiris: https://github.com/rbeyer/kalasiris
.. _ISIS: https://isis.astrogeology.usgs.gov/
.. _ISIS Interaction: https://kalasiris.readthedocs.io/en/latest/usage.html#isis-interaction


From sources
------------

The sources for hiproc can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/rbeyer/hiproc

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/rbeyer/hiproc/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/rbeyer/hiproc
.. _tarball: https://github.com/rbeyer/hiproc/tarball/master
