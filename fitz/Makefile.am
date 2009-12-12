noinst_LTLIBRARIES = libfitz.la

INCLUDES = @FREETYPE_CFLAGS@ -I$(top_srcdir)/mupdf

if HAVE_JBIG2DEC
jbig2dec_sources = filt_jbig2d.c
endif

if HAVE_OPENJPEG
openjpeg_sources = filt_jpxd.c
endif

fitz_includedir = $(includedir)
fitz_include_HEADERS =	\
	fitz.h	\
	fitz_base.h	\
	fitz_stream.h	\
	fitz_tree.h	\
	fitz_draw.h

libfitz_la_SOURCES =	\
	$(jbig2dec_sources)	\
	$(openjpeg_sources)	\
	base_cpudep.c	\
	base_error.c	\
	base_hash.c	\
	base_matrix.c	\
	base_memory.c	\
	base_rect.c	\
	base_string.c	\
	base_unicode.c	\
	\
	util_getopt.c	\
	\
	crypt_aes.c	\
	crypt_arc4.c	\
	crypt_crc32.c	\
	crypt_md5.c	\
	\
	obj_array.c	\
	obj_dict.c	\
	obj_parse.c	\
	obj_print.c	\
	obj_simple.c	\
	\
	stm_buffer.c	\
	stm_filter.c	\
	stm_open.c	\
	stm_read.c	\
	stm_misc.c	\
	\
	filt_pipeline.c	\
	filt_basic.c	\
	filt_arc4.c	\
	filt_aesd.c	\
	\
	filt_dctd.c	\
	filt_faxd.c	\
	filt_faxdtab.c	\
	filt_flate.c	\
	filt_lzwd.c	\
	filt_predict.c	\
	\
	node_toxml.c	\
	node_misc1.c	\
	node_misc2.c	\
	node_path.c	\
	node_text.c	\
	node_tree.c	\
	\
	res_colorspace.c	\
	res_font.c	\
	res_image.c	\
	res_shade.c
