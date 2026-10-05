.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in CVODES and IDAS where :c:func:`CVodeB` and :c:func:`IDASolveB`
could write past the end of the stored interpolation data when recomputing the
forward solution from a check point took more steps than the first pass.
:c:func:`CVodeB` and :c:func:`IDASolveB` now return ``CV_FWD_FAIL`` and
``IDA_FWD_FAIL``, respectively, with an error message (`Issue #49
<https://github.com/llnl/sundials/issues/49>`__).

**Deprecation Notices**
