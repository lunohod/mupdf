#include "fitz.h"

fz_shade *
fz_keepshade(fz_shade *shade)
{
	shade->refs ++;
	return shade;
}

void
fz_dropshade(fz_shade *shade)
{
	if (shade && --shade->refs == 0)
	{
		if (shade->cs)
			fz_dropcolorspace(shade->cs);
		fz_free(shade->mesh);
		fz_free(shade);
	}
}

fz_rect
fz_boundshade(fz_shade *shade, fz_matrix ctm)
{
	ctm = fz_concat(shade->matrix, ctm);
	return fz_transformrect(ctm, shade->bbox);
}

void
fz_debugshade(fz_shade *shade)
{
	int i, j, n;
	float *vert;

	printf("shade {\n");

	printf("  bbox [%g %g %g %g]\n",
		shade->bbox.x0, shade->bbox.y0,
		shade->bbox.x1, shade->bbox.y1);

	printf("  colorspace %s\n", shade->cs->name);

	printf("  matrix [%g %g %g %g %g %g]\n",
			shade->matrix.a, shade->matrix.b, shade->matrix.c,
			shade->matrix.d, shade->matrix.e, shade->matrix.f);

	if (shade->usebackground)
	{
		printf("  background [");
		for (i = 0; i < shade->cs->n; i++)
			printf("%s%g", i == 0 ? "" : " ", shade->background[i]);
		printf("]\n");
	}

	if (shade->usefunction)
	{
		printf("  function\n");
		n = 3;
	}
	else
		n = 2 + shade->cs->n;

	vert = shade->mesh;
	i = 0;
	while (i < shade->meshlen)
	{
		printf("  (%g, %g): ", vert[0], vert[1]);

		for (j = 2; j < n; j++)
			printf("%s%g", j == 2 ? "" : " ", vert[j]);
		printf("\n");

		vert += n;
		i++;
	}

	printf("}\n");
}


