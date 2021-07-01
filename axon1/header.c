#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>

/* FORWARD DECLARATIONS */
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  FILE *f;
  char tc;
  int w, h;
  int m;

  if (argc != 2)
    Error("Usage: header image.pgm\n");
  f = fopen(argv[1], "r");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open image file %s\n", argv[1]);
      exit(1);
    }
  if (!ReadHeader(f, &tc, &w, &h, &m))
    {
      fclose(f);
      fprintf(stderr, "Invalid image file header: %s\n", argv[1]);
      exit(1);
    }
  if (tc == '6')
    printf("RGB pgm\n");
  else if (tc == '5')
    printf("Grayscale pgm\n");
  else if (tc == '4')
    printf("Bilevel pgm\n");
  else
    printf("Unrecognized file type\n");
  printf("Width = %d\n", w);
  printf("Height = %d\n", h);
  if (tc != '4')
    {
      if (m == 255)
	printf("Bits = 8\n");
      else if (m == 65535)
	printf("Bits = 16\n");
      else
	printf("Unrecognized max value\n");
    }
  return(0);
}

int
ReadHeader (FILE *f, char *tc, int *w, int *h, int *m)
{
  int c;
  int v;

  c = fgetc(f);
  while (c == '#')
    {
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    }
  if (c != 'P')
    return(0);
  c = fgetc(f);
  if (c != '4' && c != '5' && c != '6')
    return(0);
  *tc = c;
  c = fgetc(f);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *w = v;
    }
  else
    return(0);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *h = v;
    }
  else
    return(0);

  if (*tc == '4')
    return(1);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *m = v;
    }
  else
    return(0);

  return(1);
}

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}
