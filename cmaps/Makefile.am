noinst_LTLIBRARIES = libcmaps.la

INCLUDES =	\
	-I$(top_srcdir)/fitz	\
	-I$(top_srcdir)/mupdf

libcmaps_la_SOURCES =	\
	cmap_tounicode.c	\
	cmap_cns.c			\
	cmap_gb.c			\
	cmap_japan.c		\
	cmap_korea.c

cmap_tounicode_files = \
	Adobe-CNS1-UCS2 Adobe-GB1-UCS2 \
	Adobe-Japan1-UCS2 Adobe-Korea1-UCS2 

cmap_cns_files = \
	Adobe-CNS1-0 Adobe-CNS1-1 Adobe-CNS1-2 Adobe-CNS1-3 \
	Adobe-CNS1-4 Adobe-CNS1-5 Adobe-CNS1-6 B5-H B5-V B5pc-H B5pc-V \
	CNS-EUC-H CNS-EUC-V CNS1-H CNS1-V CNS2-H CNS2-V ETen-B5-H \
	ETen-B5-V ETenms-B5-H ETenms-B5-V ETHK-B5-H ETHK-B5-V \
	HKdla-B5-H HKdla-B5-V HKdlb-B5-H HKdlb-B5-V HKgccs-B5-H \
	HKgccs-B5-V HKm314-B5-H HKm314-B5-V HKm471-B5-H HKm471-B5-V \
	HKscs-B5-H HKscs-B5-V UniCNS-UCS2-H UniCNS-UCS2-V 

cmap_gb_files = \
	Adobe-GB1-0 Adobe-GB1-1 Adobe-GB1-2 Adobe-GB1-3 Adobe-GB1-4 \
	Adobe-GB1-5 GB-EUC-H GB-EUC-V GB-H GB-V GBK-EUC-H GBK-EUC-V \
	GBK2K-H GBK2K-V GBKp-EUC-H GBKp-EUC-V GBpc-EUC-H GBpc-EUC-V \
	GBT-EUC-H GBT-EUC-V GBT-H GBT-V GBTpc-EUC-H GBTpc-EUC-V \
	UniGB-UCS2-H UniGB-UCS2-V 

cmap_japan_files = \
	78-EUC-H 78-EUC-V 78-H 78-RKSJ-H 78-RKSJ-V 78-V 78ms-RKSJ-H \
	78ms-RKSJ-V 83pv-RKSJ-H 90ms-RKSJ-H 90ms-RKSJ-V 90msp-RKSJ-H \
	90msp-RKSJ-V 90pv-RKSJ-H 90pv-RKSJ-V Add-H Add-RKSJ-H \
	Add-RKSJ-V Add-V Adobe-Japan1-0 Adobe-Japan1-1 Adobe-Japan1-2 \
	Adobe-Japan1-3 Adobe-Japan1-4 Adobe-Japan1-5 Adobe-Japan1-6 \
	EUC-H EUC-V Ext-H Ext-RKSJ-H Ext-RKSJ-V Ext-V H Hankaku \
	Hiragana Katakana NWP-H NWP-V RKSJ-H RKSJ-V Roman \
	UniJIS-UCS2-H UniJIS-UCS2-HW-H UniJIS-UCS2-HW-V UniJIS-UCS2-V \
	UniJISPro-UCS2-HW-V UniJISPro-UCS2-V V WP-Symbol \
	Adobe-Japan2-0 Hojo-EUC-H Hojo-EUC-V Hojo-H Hojo-V \
	UniHojo-UCS2-H UniHojo-UCS2-V 

cmap_korea_files =  \
	Adobe-Korea1-0 Adobe-Korea1-1 Adobe-Korea1-2 KSC-EUC-H \
	KSC-EUC-V KSC-H KSC-Johab-H KSC-Johab-V KSC-V KSCms-UHC-H \
	KSCms-UHC-HW-H KSCms-UHC-HW-V KSCms-UHC-V KSCpc-EUC-H \
	KSCpc-EUC-V UniKS-UCS2-H UniKS-UCS2-V 

cmap_tounicode.c: $(cmap_tounicode_files)
	cmapdump $@ $+

cmap_cns.c: $(cmap_cns_files)
	cmapdump $@ $+

cmap_gb.c: $(cmap_gb_files)
	cmapdump $@ $+

cmap_japan.c: $(cmap_japan_files)
	cmapdump $@ $+

cmap_korea.c: $(cmap_korea_files)
	cmapdump $@ $+
