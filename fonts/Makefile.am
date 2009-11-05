noinst_LTLIBRARIES = libfonts.la

if BUILD_CJK
cjk_sources = font_cjk.c
endif

INCLUDES =	\
	-I$(top_srcdir)/fitz	\
	-I$(top_srcdir)/mupdf

libfonts_la_SOURCES =	\
	$(cjk_sources)  \
	font_misc.c	\
	font_mono.c	\
	font_serif.c	\
	font_sans.c

font_misc.c: Dingbats.cff StandardSymL.cff URWChanceryL-MediItal.cff
	fontdump $@ $+

font_mono.c: NimbusMonL-Regu.cff NimbusMonL-ReguObli.cff NimbusMonL-Bold.cff NimbusMonL-BoldObli.cff
	fontdump $@ $+

font_serif.c: NimbusRomNo9L-Regu.cff NimbusRomNo9L-ReguItal.cff NimbusRomNo9L-Medi.cff NimbusRomNo9L-MediItal.cff
	fontdump $@ $+

font_sans.c: NimbusSanL-Bold.cff NimbusSanL-BoldItal.cff NimbusSanL-Regu.cff NimbusSanL-ReguItal.cff
	fontdump $@ $+

font_cjk.c: droid/DroidSansFallback.ttf
	fontdump $@ $+
	
