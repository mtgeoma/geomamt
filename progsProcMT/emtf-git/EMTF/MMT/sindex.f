c________________________________________________
c
        function sindex(i,j)
 
        integer sindex
        i1=min(i,j)
        j1=max(i,j)
        sindex=((j1-1)*j1)/2+i1
        return
        end
