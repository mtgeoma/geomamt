c   nfilm: maximum number of frequency that will be read from calibration file

c   from fitpack.f:
c   sigma contains the tension factor. this value indicates
c   the curviness desired. if abs(sigma) is nearly zero
c   (e.g. .001) the resulting curve is approximately a
c   cubic spline. if abs(sigma) is large (e.g. 50.) the
c   resulting curve is nearly a polygonal line. if sigma
c   equals zero a cubic spline results.  a standard value
c   for sigma is approximately 1. in absolute value.

      integer nfilm
      real sigma
      parameter (nfilm=100)
      parameter (sigma=0.)
