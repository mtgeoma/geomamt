c****************************************************************
      real function convz2p( z, per )
c
c   convert Z to phase
c
c****************************************************************
      complex z
      real per
 
      convz2p = atan2d( aimag(z), real(z) )
 
      return
      end
