How to install spglib C-API
============================

Download
---------

The source code is downloaded at
https://github.com/atztogo/spglib/archive/master.zip or you can
git-clone the spglib repository by

::

   % git clone https://github.com/atztogo/spglib.git

Install
--------

Compiling using cmake
^^^^^^^^^^^^^^^^^^^^^^

1. After expanding source code, go into the source code directory::

     % tar xvfz spglib-1.9.8.tar.gz
     % cd spglib-1.9.8

2. Build and install in ``_build`` directory by

   ::

     % mkdir _build && _build
     % cmake -DCMAKE_INSTALL_PREFIX="" ..
     % make
     % make DESTDIR=/some/where install

   The libraries are installed at ``/some/where/lib`` and the
   header file is installed at ``/some/where/include``.


Compiling using configure script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. The configure script is prepared using
   autotools and libtool as follows::

     % aclocal
     % autoheader
     % libtoolize # or glibtoolize with macport etc
     % touch INSTALL NEWS README AUTHORS
     % automake -acf
     % autoconf


2. Run configure script::

     % tar xvfz spglib-1.9.8.tar.gz
     % cd spglib-1.9.8
     % ./configure --prefix=INSTALLATION_LOCATION
     % make
     % make install

3. The libraries are installed at ``INSTALLATION_LOCATION/lib`` or found in
   ``src/.libs`` if you don't run ``make install``.

Usage
------

1. Include ``spglib.h``
2. Link ``libsymspg.a`` or ``libsymspg.so``
3. A compilation example is shown in  `example/README
   <https://github.com/atztogo/spglib/blob/master/example/README>`_.

Example
--------

A few examples are found in `example
<https://github.com/atztogo/spglib/tree/master/example>`_ directory.
