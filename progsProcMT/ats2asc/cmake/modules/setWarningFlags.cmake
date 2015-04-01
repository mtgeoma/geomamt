include(CheckCXXCompilerFlag)

function(setWarningFlags warningFlags_ref)
  set(${warningFlags_ref} "" PARENT_SCOPE)

  #Se CMAKE_CXX_FLAGS contiver "-Werror", os testes não funcionam.
  string(REGEX REPLACE "-Werror" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

  set(possibleWarnings
    -Wabi
    -Wconversion-null
    -Wctor-dtor-privacy
    #-Wnoexcept (-Wno)
    #-Wnon-virtual-dtor (-Weffc++)
    #-Wreorder (-Wall)
    #-Weffc++ (cmake falha)
    -Wstrict-null-sentinel
    #-Wno-non-template-friend (-Wno)
    -Wold-style-cast
    -Woverloaded-virtual
    -Wno-pmf-conversions
    -Wsign-promo
    #-Wassign-intercept (Objective-C and Objective-C++ only)
    #-Wno-protocol (Objective-C and Objective-C++ only)
    #-Wselector (Objective-C and Objective-C++ only)
    #-Wstrict-selector-match (Objective-C and Objective-C++ only)
    #-Wundeclared-selector (Objective-C and Objective-C++ only)
    #-fsyntax-only (Check the code for syntax errors only)
    #-fmax-errors=@var{n} (Limits the maximum number of error messages to n)
    -pedantic
    #-pedantic-errors (transform -pedantic warnings into erros)
    #-w (Inhibit all warning messages)
    -Wextra
    -Wall
    #-Waddress (-Wall)
    #-Waggregate-return
    #-Warray-bounds (-Wall && -O2)
    #-Wno-attributes (-Wno)
    #-Wno-builtin-macro-redefined (-Wno)
    #-Wc++-compat (C and Objective-C only) (falhou no cmake)
    #-Wc++0x-compat (-Wall)
    -Wcast-align
    -Wcast-qual
    #-Wchar-subscripts (-Wall)
    #-Wclobbered (-Wextra)
    #-Wcomment (-Wall)
    -Wconversion
    -Wcoverage-mismatch
    #-Wno-cpp (-Wno)
    #-Wno-deprecated (-Wno)
    #-Wno-deprecated-declarations (-Wno)
    -Wdisabled-optimization
    #-Wno-div-by-zero (-Wno)
    -Wdouble-promotion
    #-Wempty-body (-Wextra)
    #-Wenum-compare (-Wall)
    #-Wno-endif-labels (-Wno)
    #-Werror (Make all warnings into errors.)
    #-Werror=* (Make the specified warning into an error)
    #-Wfatal-errors (Abort compilation on the first error)
    -Wfloat-equal
    #-Wformat (-Wall) (-Wformat=2)
    -Wformat=2
    #-Wno-format-contains-nul (-Wno)
    #-Wno-format-extra-args (-Wno)
    #-Wformat-nonliteral (-Wformat=2)
    #-Wformat-security (-Wformat=2)
    #-Wformat-y2k (-Wformat=2)
    #-Wframe-larger-than=0
    #-Wno-free-nonheap-object (-Wno)
    #-Wjump-misses-init (C, Objective-C only)
    #-Wignored-qualifiers (-Wextra)
    #-Wimplicit (-Wall)
    #-Wimplicit-function-declaration (-Wall)
    #-Wimplicit-int (-Wall)
    -Winit-self
    -Winline
    #-Wmaybe-uninitialized (-Wall)
    #-Wno-int-to-pointer-cast (-Wno)
    #-Wno-invalid-offsetof (-Wno)
    -Winvalid-pch
    #-Wlarger-than=0
    -Wunsafe-loop-optimizations
    -Wlogical-op
    #-Wlong-long (-pedantic)
    #-Wmain (-Wall)
    #-Wmaybe-uninitialized (-Wall)
    #-Wmissing-braces (-Wall)
    #-Wmissing-field-initializers (-Wextra)
    -Wmissing-format-attribute
    -Wmissing-include-dirs
    #-Wno-mudflap (-Wno)
    #-Wno-multichar (-Wno)
    #-Wnonnull (-Wall)
    #-Wno-overflow (-Wno)
    #-Woverlength-strings (-pedantic)
    -Wpacked
    -Wpacked-bitfield-compat
    -Wpadded
    #-Wparentheses (-Wall)
    -Wpedantic-ms-format
    #-Wno-pedantic-ms-format (-Wno)
    #-Wpointer-arith (-pedantic)
    #-Wno-pointer-to-int-cast (-Wno)
    -Wredundant-decls
    #-Wreturn-type (-Wall)
    #-Wsequence-point (-Wall)
    -Wshadow
    #-Wsign-compare (-Wall) (-Wextra)
    -Wsign-conversion
    -Wstack-protector
    -Wstack-usage=0
    #-Wstrict-aliasing (-Wall)
    #-Wstrict-aliasing=n (-Wstrict-aliasing)
    #-Wstrict-overflow (-Wstrict-overflow=)
    -Wstrict-overflow=5
    -Wsuggest-attribute=pure
    -Wsuggest-attribute=const
    -Wsuggest-attribute=noreturn
    #-Wswitch (-Wall)
    -Wswitch-default
    -Wswitch-enum
    -Wsync-nand
    #-Wsystem-headers (Warning messages found in system header files)
    -Wtrampolines
    #-Wtrigraphs (-Wall)
    #-Wtype-limits (-Wextra)
    -Wundef
    #-Wuninitialized (-Wall) (-Wextra)
    #-Wunknown-pragmas (-Wall)
    #-Wno-pragmas (-Wno)
    #-Wunsuffixed-float-constants (C and Objective-C only)
    -Wunused
    #-Wunused-function (-Wall) || (-Wunused)
    #-Wunused-label (-Wall) || (-Wunused)
    #-Wunused-local-typedefs (-Wunused)
    #-Wunused-parameter (-Wunused) || (-Wextra && -Wall)
    #-Wno-unused-result (-Wno)
    #-Wunused-value (-Wall) || (-Wunused)
    #-Wunused-variable (-Wall) (-Wunused)
    #-Wunused-but-set-parameter (-Wunused) || (-Wextra && -Wall)
    #-Wunused-but-set-variable (-Wunused)
    -Wvariadic-macros
    -Wvector-operation-performance
    -Wvla
    #-Wvolatile-register-var (-Wall)
    -Wwrite-strings
    )

  foreach(flag IN LISTS possibleWarnings)
    check_cxx_compiler_flag(${flag} Warning${flag})
    if(${Warning${flag}})
      set(warningFlags "${warningFlags} ${flag}")
    endif()
  endforeach()

  #-Weffc++ não funciona no loop acima
  check_cxx_compiler_flag(-Weffc++ Warning-Weffcpp)
  if(${Warning-Weffcpp})
    set(warningFlags "${warningFlags} -Weffc++")
  endif()

  set(${warningFlags_ref} "${warningFlags}" PARENT_SCOPE)
endfunction()

setWarningFlags(cxxWarningFlags)
