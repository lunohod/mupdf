lib_LTLIBRARIES = libmupdf.la

if BUILD_CJK
else
AM_CFLAGS = -DNOCJK
endif

INCLUDES = @FREETYPE_CFLAGS@ @JASPER_CFLAGS@ @JBIG2DEC_CFLAGS@ @OPENJPEG_CFLAGS@ -I$(top_srcdir)/fitz

mupdf_includedir = $(includedir)
mupdf_include_HEADERS = mupdf.h

libmupdf_la_LIBADD =	\
	$(top_builddir)/fitz/libfitz.la	\
	$(top_builddir)/fitzdraw/libfitzdraw.la	\
	$(top_builddir)/fonts/libfonts.la	\
	$(top_builddir)/cmaps/libcmaps.la	\
	@FREETYPE_LIBS@	\
	@JASPER_LIBS@ \
	@JBIG2DEC_LIBS@ \
	@OPENJPEG_LIBS@ \
	@JPEG_LIBS@ \
	@Z_LIBS@ \
	-lm
	
libmupdf_la_LDFLAGS = -no-undefined -version-info 1:0:0 -Wl,--as-needed
libmupdf_la_SOURCES =	\
	pdf_crypt.c	\
	pdf_debug.c	\
	pdf_lex.c	\
	pdf_nametree.c	\
	pdf_open.c	\
	pdf_parse.c	\
	pdf_repair.c	\
	pdf_stream.c	\
	pdf_xref.c	\
	\
	pdf_annot.c	\
	pdf_outline.c	\
	\
	pdf_cmap.c	\
	pdf_cmap_parse.c	\
	pdf_cmap_load.c	\
	pdf_cmap_table.c	\
	\
	pdf_fontagl.c	\
	pdf_fontenc.c	\
	pdf_unicode.c	\
	pdf_font.c	\
	pdf_type3.c	\
	pdf_fontmtx.c	\
	pdf_fontfile.c	\
	\
	pdf_function.c	\
	pdf_colorspace1.c	\
	pdf_colorspace2.c	\
	pdf_image.c	\
	pdf_pattern.c	\
	pdf_shade.c	\
	pdf_shade1.c	\
	pdf_shade4.c	\
	pdf_xobject.c	\
	\
	pdf_build.c	\
	pdf_interpret.c	\
	pdf_page.c	\
	pdf_pagetree.c	\
	pdf_store.c


