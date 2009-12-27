#include "fitz_base.h"
#include "fitz_tree.h"
#include "fitz_draw.h"

typedef unsigned char byte;

static fz_error gettile(fz_lazytile *l, int y0)
{
	fz_error error;
	fz_image *image = l->image;
	int dx = l->dx;
	int dy = l->dy;
	fz_pixmap *tile;
	fz_pixmap *temp;
	int tileheight = l->tileheight;

	if(y0 < 0)
		y0 = 0;

	int y = y0 * dy;

	if (y + tileheight > image->h)
		tileheight = image->h - y;

	if (tileheight < 0)
		return fz_okay;

	if (l->tile)
	{
		fz_droppixmap(l->tile);
		l->tile = NULL;
	}

	error = fz_newpixmap(&tile, 0, y, image->w, tileheight, image->n + 1);
	if (error)
		return error;

	error = image->loadtile(image, tile);
	if (error)
		goto cleanup;

	if (dx > 1 || dy > 1)
	{
		error = fz_scalepixmap(&temp, tile, dx, dy);
		if (error)
			goto cleanup;

		fz_droppixmap(tile);
		tile = temp;
	}

	if (image->cs && image->cs != l->target_cs)
	{
		error = fz_newpixmap(&temp, tile->x, tile->y, tile->w, tile->h, l->target_cs->n + 1);
		if (error)
			goto cleanup;

		fz_convertpixmap(image->cs, tile, l->target_cs, temp);
		fz_droppixmap(tile);
		tile = temp;
	}

	tile->y = y0;
	l->tile = tile;

	return fz_okay;

cleanup:
	if (tile)
		fz_droppixmap(tile);

	return error;
}

static inline void checktile(fz_lazytile *l, int y)
{
	if (!l->tile)
	{
		gettile(l, y);
		return;
	}

	int y0 = l->tile->y;

	while (y0 > y)
		y0 -= (l->tileheight / l->dy) + 1;
	while (y0 + l->tile->h <= 1 + y)
		y0 += (l->tileheight / l->dy) - 1;

	if (y0 < 0)
		y0 = 0;

	if (l->tile->y != y0)
		gettile(l, y0);
}

#define lerp(a,b,t) (a + (((b - a) * t) >> 16))

static inline byte getcomp(fz_lazytile *s, int w, int h, int u, int v, int n, int k)
{
	if (u < 0) u = 0;
	if (v < 0) v = 0;
	if (u >= w) u = w - 1;
	if (v >= h) v = h - 1;

	return s->tile->samples[(w * (v - s->tile->y) + u) * n + k];
}

static inline int samplecomp(fz_lazytile *s, int w, int h, int u, int v, int n, int k)
{
	int ui = u >> 16;
	int vi = v >> 16;
	int ud = u & 0xFFFF;
	int vd = v & 0xFFFF;
	int a = getcomp(s, w, h, ui, vi, n, k);
	int b = getcomp(s, w, h, ui+1, vi, n, k);
	int c = getcomp(s, w, h, ui, vi+1, n, k);
	int d = getcomp(s, w, h, ui+1, vi+1, n, k);
	int ab = lerp(a, b, ud);
	int cd = lerp(c, d, ud);
	return lerp(ab, cd, vd);
}

static inline byte getmask(fz_lazytile *s, int w, int h, int u, int v)
{
	if (u < 0) u = 0;
	if (v < 0) v = 0;
	if (u >= w) u = w - 1;
	if (v >= h) v = h - 1;

	return s->tile->samples[w * (v - s->tile->y) + u];
}

static inline int samplemask(fz_lazytile *s, int w, int h, int u, int v)
{
	int ui = u >> 16;
	int vi = v >> 16;
	int ud = u & 0xFFFF;
	int vd = v & 0xFFFF;
	int a = getmask(s, w, h, ui, vi);
	int b = getmask(s, w, h, ui+1, vi);
	int c = getmask(s, w, h, ui, vi+1);
	int d = getmask(s, w, h, ui+1, vi+1);
	int ab = lerp(a, b, ud);
	int cd = lerp(c, d, ud);
	return lerp(ab, cd, vd);
}

static inline void lerpargb(byte *dst, byte *a, byte *b, int t)
{
	dst[0] = lerp(a[0], b[0], t);
	dst[1] = lerp(a[1], b[1], t);
	dst[2] = lerp(a[2], b[2], t);
	dst[3] = lerp(a[3], b[3], t);
}

static inline byte *getargb(fz_lazytile *s, int w, int h, int u, int v)
{
	if (u < 0) u = 0;
	if (v < 0) v = 0;
	if (u >= w) u = w - 1;
	if (v >= h) v = h - 1;

	return s->tile->samples + ((w * (v - s->tile->y) + u) << 2);
}

static inline void sampleargb(fz_lazytile *s, int w, int h, int u, int v, byte *abcd)
{
	byte ab[4];
	byte cd[4];
	int ui = u >> 16;
	int vi = v >> 16;
	int ud = u & 0xFFFF;
	int vd = v & 0xFFFF;
	byte *a = getargb(s, w, h, ui, vi);
	byte *b = getargb(s, w, h, ui+1, vi);
	byte *c = getargb(s, w, h, ui, vi+1);
	byte *d = getargb(s, w, h, ui+1, vi+1);
	lerpargb(ab, a, b, ud);
	lerpargb(cd, c, d, ud);
	lerpargb(abcd, ab, cd, vd);
}

static void img_ncn(FZ_LPSRC, int srcn, FZ_PDST, FZ_PCTM)
{
	int k;
	while (h--)
	{
		byte *dstp = dst0;
		int u = u0;
		int v = v0;
		int w = w0;

		checktile(src, v>>16);

		while (w--)
		{
			for (k = 0; k < srcn; k++)
			{
				dstp[k] = samplecomp(src, srcw, srch, u, v, srcn, k);
				dstp += srcn;
				u += fa;
				v += fb;
			}
		}
		dst0 += dstw;
		u0 += fc;
		v0 += fd;
	}
}

static void img_1c1(FZ_LPSRC, FZ_PDST, FZ_PCTM)
{
	while (h--)
	{
		byte *dstp = dst0;
		int u = u0;
		int v = v0;
		int w = w0;

		checktile(src, v>>16);

		while (w--)
		{
			dstp[0] = samplemask(src, srcw, srch, u, v);
			dstp ++;
			u += fa;
			v += fb;
		}
		dst0 += dstw;
		u0 += fc;
		v0 += fd;
	}
}

static void img_4c4(FZ_LPSRC, FZ_PDST, FZ_PCTM)
{
	while (h--)
	{
		byte *dstp = dst0;
		int u = u0;
		int v = v0;
		int w = w0;

		checktile(src, v>>16);

		while (w--)
		{
			sampleargb(src, srcw, srch, u, v, dstp);
			dstp += 4;
			u += fa;
			v += fb;
		}
		dst0 += dstw;
		u0 += fc;
		v0 += fd;
	}
}

static void img_1o1(FZ_LPSRC, FZ_PDST, FZ_PCTM)
{
	byte srca;
	while (h--)
	{
		byte *dstp = dst0;
		int u = u0;
		int v = v0;
		int w = w0;

		checktile(src, v>>16);

		while (w--)
		{
			srca = samplemask(src, srcw, srch, u, v);
			dstp[0] = srca + fz_mul255(dstp[0], 255 - srca);
			dstp ++;
			u += fa;
			v += fb;
		}
		dst0 += dstw;
		u0 += fc;
		v0 += fd;
	}
}

static void img_4o4(FZ_LPSRC, FZ_PDST, FZ_PCTM)
{
	byte argb[4];
	byte ssa;
	while (h--)
	{
		byte *dstp = dst0;
		int u = u0;
		int v = v0;
		int w = w0;

		checktile(src, v>>16);

		while (w--)
		{
			sampleargb(src, srcw, srch, u, v, argb);
			ssa = 255 - argb[0];
			dstp[0] = argb[0] + fz_mul255(dstp[0], ssa);
			dstp[1] = argb[1] + fz_mul255(dstp[1], ssa);
			dstp[2] = argb[2] + fz_mul255(dstp[2], ssa);
			dstp[3] = argb[3] + fz_mul255(dstp[3], ssa);
			dstp += 4;
			u += fa;
			v += fb;
		}
		dst0 += dstw;
		u0 += fc;
		v0 += fd;
	}
}

static void img_w4i1o4(byte *argb, FZ_LPSRC, FZ_PDST, FZ_PCTM)
{
	byte alpha = argb[0];
	byte r = argb[4];
	byte g = argb[5];
	byte b = argb[6];
	byte cov;
	byte ca;
	while (h--)
	{
		byte *dstp = dst0;
		int u = u0;
		int v = v0;
		int w = w0;

		checktile(src, v>>16);

		while (w--)
		{
			cov = samplemask(src, srcw, srch, u, v);
			ca = fz_mul255(cov, alpha);
			dstp[0] = ca + fz_mul255(dstp[0], 255 - ca);
			dstp[1] = fz_mul255((short)r - dstp[1], ca) + dstp[1];
			dstp[2] = fz_mul255((short)g - dstp[2], ca) + dstp[2];
			dstp[3] = fz_mul255((short)b - dstp[3], ca) + dstp[3];
			dstp += 4;
			u += fa;
			v += fb;
		}
		dst0 += dstw;
		u0 += fc;
		v0 += fd;
	}
}

void (*fz_img_ncn_lazy)(FZ_LPSRC, int sn, FZ_PDST, FZ_PCTM) = img_ncn;
void (*fz_img_1c1_lazy)(FZ_LPSRC, FZ_PDST, FZ_PCTM) = img_1c1;
void (*fz_img_4c4_lazy)(FZ_LPSRC, FZ_PDST, FZ_PCTM) = img_4c4;
void (*fz_img_1o1_lazy)(FZ_LPSRC, FZ_PDST, FZ_PCTM) = img_1o1;
void (*fz_img_4o4_lazy)(FZ_LPSRC, FZ_PDST, FZ_PCTM) = img_4o4;
void (*fz_img_w4i1o4_lazy)(byte*,FZ_LPSRC,FZ_PDST,FZ_PCTM) = img_w4i1o4;
