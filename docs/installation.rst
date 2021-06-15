.. highlight:: shell

============
Installation
============


.. Stable release
   --------------

To install hiproc, you must first do some messing around.  Someday, all of
this dependency madness might be smoothed out.  Of course, it could always be
worse.

Firstly, this package depends on `kalasiris`_ which is just a wrapper
around `ISIS`_, which has some complicated interactions with other
environments.  You will need to install ISIS separately in another
conda environment (we'll call it *isis* for the purpose of this
discussion), and then review the strategies detailed in the kalasiris
page on `ISIS Interaction`_ to ensure that you can get everything
working.

At any rate, create the *isis* conda environment, run the init script,
install the data directories that you need, etc.  Odds are good that you
probably already have an *isis* environment, so maybe you're already good.

Since you've already got conda installed, the remainder of these instructions
are going to use conda (you could also use pip or something else).

Make a new conda (or other kind of virtual) environment, let's call it
*hiproc*.

Now, you're finally ready to install hiproc (make sure the conda-forge channel is
available):

    .. code-block:: console

        $ conda install hiproc

You could also install `from sources`_ at this point, instead.

So now we have hiproc in the *hiproc* environment, but a separate *isis* environment.

The easiest thing to do is exit out of all conda envrionments, and then:

    .. code-block:: console

        $ conda activate isis
        $ conda activate --stack hiproc

Then you should be ready to run hiproc, and it will be able to access your ISIS install.

If you are going to try and use hiproc with the Ames Stereo Pipeline (ASP), you'll do something
a little different.  We'll assume you have followed those instructions and have an *asp* environment
(which installs its own ISIS), and a *hiproc* environment.  In this case, you wouldn't need a
separate *isis* environment, but you might have it anyway.  To get things working together,
exit out of any conda environments, and do this:

    .. code-block:: console

        $ conda activate hiproc
        $ conda activate --stack asp

And now, you should have access to ISIS and ASP programs via the *asp* environment, but also
hiproc via the *hiproc* environment.  Someday, when the ISIS dependencies aren't so strict,
we'll be able to install all of this in a single conda environment, and we won't have to
play these games.

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


On a system with the Perl versions
----------------------------------
In the very rare case that you are working on a system where the
Perl versions of these programs are in the default path (this likely means
that you are working for HiRISE or on a HiRISE Team system), installation
of this Python package may shadow those paths.  For example, running
"HiCal" may run this Python version, and *not* the Perl version you might
have been expecting.

You can alleviate this by always using absolute paths (tedious), or you can
do some tweaking so that the Python versions are either not installed on the path
or are installed with different names that you control.

You will have to follow the above instructions, and then use the `from sources`_
method to get a copy of the GitHub repo.

Then, before "installing" it with the ``python setup.py install``,
you will need to edit ``setup.py``.

In ``setup.py`` there is a section where a dict is defined::

    entry_points={
        'console_scripts': [
            'EDR_Stats=hiproc.EDR_Stats:main',
            'HiBeautify=hiproc.HiBeautify:main',

Here, you should delete the lines of the programs that you don't want installed.  Or
you could change them to give them different names, so ``'EDR_Stats=hiproc.EDR_Stats:main'``
could become ``'hiproc-EDR_Stats=hiproc.EDR_Stats:main'`` and then after "installation" you'd
have a program named  ``hiproc-EDR_Stats`` installed in your path, and your pre-exsiting
Perl ``EDR_Stats`` would still be accessible.

After you edit the ``setup.py`` file, then you can ``python setup.py install`` or
``pip install --no-deps -e .`` whatever you prefer to get it "installed."

.. _Github repo: https://github.com/rbeyer/hiproc
.. _tarball: https://github.com/rbeyer/hiproc/tarball/master
