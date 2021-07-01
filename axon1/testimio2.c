#include <stdio.h>
#include <stdlib.h>
#include "imio.h"

int
main ()
{
  unsigned char *img;
  char msg[4096];
  int w, h;

  if (!ReadImage("im", &img, &w, &h, -1, -1, -1, -1, msg))
    {
      fprintf(stderr, "ReadImage: %s\n", msg);
      exit(1);
    }
  if (!WriteImage("im.pgm", img, w, h, UncompressedImage, msg))
    {
      fprintf(stderr, "WriteImage: %s\n", msg);
      exit(1);
    }
    
  return(0);
}
