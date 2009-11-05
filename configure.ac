AC_PREREQ(2.57)
AC_INIT(mupdf, 0.4, http://bugs.ghostscript.com/)

AM_INIT_AUTOMAKE(foreign)

AC_PROG_CC_C99
AC_PROG_LIBTOOL

PKG_CHECK_MODULES(FREETYPE, freetype2)

AC_CHECK_LIB(jasper, jas_init, [JASPER_LIBS=-ljasper JASPER_CFLAGS=-DHAVE_JASPER], [], -ljpeg)
AC_SUBST(JASPER_LIBS)
AC_SUBST(JASPER_CFLAGS)
AM_CONDITIONAL(HAVE_JASPER, test "x$JASPER_LIBS" != x)

AC_CHECK_LIB(jbig2dec, jbig2_ctx_new, [JBIG2DEC_LIBS=-ljbig2dec JBIG2DEC_CFLAGS=-DHAVE_JBIG2DEC], [], [])
AC_SUBST(JBIG2DEC_LIBS)
AC_SUBST(JBIG2DEC_CFLAGS)
AM_CONDITIONAL(HAVE_JBIG2DEC, test "x$JBIG2DEC_LIBS" != x)

AC_CHECK_LIB(openjpeg, opj_cio_open, [OPENJPEG_LIBS=-lopenjpeg OPENJPEG_CFLAGS=-DHAVE_OPENJPEG], [], [])
AC_SUBST(OPENJPEG_LIBS)
AC_SUBST(OPENJPEG_CFLAGS)
AM_CONDITIONAL(HAVE_OPENJPEG, test "x$OPENJPEG_LIBS" != x)

AC_CHECK_LIB(jpeg, jpeg_read_header, [JPEG_LIBS=-ljpeg], AC_MSG_ERROR([Libjpeg library is not available]), [])
AC_SUBST(JPEG_LIBS)

AC_CHECK_LIB(z, zlibVersion, [Z_LIBS=-lz], AC_MSG_ERROR([zlib library is not available]), [])
AC_SUBST(Z_LIBS)

AC_ARG_ENABLE(cjk, 
	AC_HELP_STRING([--enable-cjk],
			[enable builtin CJK fonts and data]),
		enable_cjk=$enableval,
		enable_cjk="no")
AM_CONDITIONAL(BUILD_CJK, test x$enable_cjk = xyes)

AC_OUTPUT([
Makefile
mupdf.pc
fitz/Makefile
fitzdraw/Makefile
cmaps/Makefile
cmaps/chinese_s/Makefile
cmaps/chinese_t/Makefile
cmaps/japanese/Makefile
cmaps/korean/Makefile
fonts/Makefile
mupdf/Makefile
])
