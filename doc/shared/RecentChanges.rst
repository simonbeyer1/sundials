.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed bugs in CVODES and IDAS where the interpolation of the forward solution
for adjoint sensitivity analysis could return a wrong value without an error or
read past the end of the stored data. With polynomial interpolation, this
happened when a check point interval held fewer steps than the interpolation
order or when integrating toward decreasing time. With either interpolation
type, it happened after an evaluation at the start of a check point interval.
:c:func:`IDAInitB`, :c:func:`IDAInitBS`, and :c:func:`IDAReInitB` now accept
``tB0`` when integrating toward decreasing time.

**Deprecation Notices**
