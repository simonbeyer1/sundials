.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in CVODES and IDAS where the polynomial interpolation of the forward
solution for adjoint sensitivity analysis could return a wrong value without an
error or read past the end of the stored data when a check point interval held
fewer steps than the interpolation order.

**Deprecation Notices**
