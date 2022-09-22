#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <njansson@kth.se> wrote this file. As long as you retain this notice you
# can do whatever you want with this stuff. If we meet some day, and you think
# this stuff is worth it, you can buy me a beer in return Niclas Jansson
# ----------------------------------------------------------------------------
#
AC_DEFUN([AX_KAHIP],[
        AC_ARG_WITH([kahip],
                    AC_HELP_STRING([--with-kahip=DIR],
                    [Compile with support for KaHIP]),
                    [
                    if test -d "$withval"; then
                       ac_kahip_path="$withval";
                       AS_IF([test -d "$ac_kahip_path/lib64"],
                             [suffix="64"],[suffix=""])
                       KAHIP_LDFLAGS="-L$ac_kahip_path/lib$suffix"
                       KAHIP_CPPFLAGS="-I$ac_kahip_path/include"
                    fi
                    ], [with_kahip=no])
                    
        if test "x${with_kahip}" != xno; then

           if test -d "$ac_kahip_path"; then
	      CPPFLAGS_SAVED="$CPPFLAGS"
              LDFLAGS_SAVED="$LDFLAGS"
	      CPPFLAGS="$KAHIP_CPPFLAGS $CPPFLAGS"
	      LDFLAGS="$KAHIP_LDFLAGS $LDFLAGS"
	      export CPPFLAGS
	      export LDFLAGS
	   fi

           _LIBS=$LIBS
           LIBS=""
	   AC_LANG_PUSH([C++])
	   AC_CHECK_HEADER([parhip_interface.h],
                           [have_parhip_interface_h=yes],
                           [have_parhip_interface_h=no])
	   if test x"${have_parhip_interface_h}" = xno; then
		if test -d "$ac_kahip_path"; then	
		   CPPFLAGS="$CPPFLAGS_SAVED"
		fi
	   fi
           AC_LANG_POP([C++])

           AC_LANG_PUSH([C])   

	   AC_CHECK_LIB(parhip_interface, ParHIPPartitionKWay,
			[have_kahip=yes;KAHIP_LIBS="-lparhip_interface"],
			[have_kahip=no],[])
	   AC_SUBST(KAHIP_LIBS)
	   if test x"${have_kahip}" = xyes; then
	      AC_DEFINE(HAVE_KAHIP,1,[Define if you have the KaHIP library.])
              LIBS="$KAHIP_LIBS $_LIBS"
	   else
	      if test -d "$ac_kahip_path"; then	
	      	 LDFLAGS="$LDFLAGS_SAVED"
	      fi
	   fi

           AC_LANG_POP([C])   
       fi
])
