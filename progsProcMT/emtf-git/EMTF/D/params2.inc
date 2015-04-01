c       include file PARAMS2.INC
c
c       This include file is included by DECIMATE.INC

c       PARAMETERS:
c
c       NDMX:   Maximum number of decimation levels
c       NCHMX:  Maximum number of channels to be fft'd
c       NFCMX:  Maximum number of filter coefficients for
c               decimation filter (should be about equal 
c               to the maximum decimation factor)
c       NFILMAX:  Maximum number of filters to correct for for one
c               data channel

        parameter       (ndmx = 6)
        include '../include/nchmx.inc'
        parameter       (nfcmx = 5)
        parameter       (nfilmax = 4)
