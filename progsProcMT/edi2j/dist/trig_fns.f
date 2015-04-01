c      routines for trig functions not found on some machines

      function sind(x)
      real x
      sind = sin(x*3.14159265/180.0)
      return
      end

      function cosd(x)
      real x
      cosd = cos(x*3.14159265/180.0)
      return
      end
                  
      function tand(x)
      real x
      tand = tan(x*3.14159265/180.0)
      return
      end

      function asind(x)
      real x
      asind = (180.0/3.14159265)*sin(x)
          return
      end

      function acosd(x)
      real x
      acosd = (180.0/3.14159265)*cos(x)
          return
      end

      function atand(x)
      real x
      atand = (180.0/3.14159265)*atan(x)
      return
      end
            
      double precision function dtand(x)
      double precision x
      dtand = dtan(x*3.14159265d0/180.d0)
      return
      end
                   
      function atan2d(y,x)
      real x,y
      atan2d = (180.0/3.14159265)*atan2(y,x)
      return
      end     
            
      double precision function datand(x)
      double precision x
      datand = (180.d0/3.14159265d0)*datan(x)
      return
      end
            
      double precision function datan2d(y,x)
      double precision x,y
      datan2d = (180.d0/3.14159265d0)*datan2(y,x)
      return
      end
