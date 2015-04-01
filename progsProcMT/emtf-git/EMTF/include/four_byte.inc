      logical l_4byte
ccc   set l_4byte = .true. if FOrtran compiler assumes that record length for
ccc   direct access file connect is specified in four byte words (instead of
ccc   bytes)
      data l_4byte/.false./
