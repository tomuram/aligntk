/* APPENDIX

Code Fragment for Generating the Two Dimensional Hilbert Mapping
*/

int rotation_table[4] = {3, 0, 0, 1};
int sense_table[4] = {-1, 1, 1, -1};

int quad_table[4][2][2]= { {{0,1},{3,2}},
			   {{1,2},{0,3}},{{2,3},{1,0}},{{3,0},{2,1}} };

{
  rotation = 0; /* Initially no rotation */
  sense = 1; /* Initially positive sense */
  num=O;
  for (k=side/2; k>O; k=k/2)
    {
      xbit = x/k; /*Get the msb of x*/
      ybit = y/k;
      x -= k*xbit; /* Take away the current m s b */
      y -= k*ybit;
      quad = quad_table[rotation][xbit][ybit]; /*Which quadrant am I in? */
      num += (sense == - 1) ? k*k*(3-quad) : k*k*quad;
      rotation += rotation_table[quad]; /* Fix rotation value for next time */
      if (rotation >= 4) rotation -= 4; /* Addition is modulo 4 */
      sense *= sense_table[quad]; /* Fix sense value for next time */
    }
}
/* The above code fragment produces in num the linear
   position of the grid point with coordinates x and y, in square
   grid with side of size side */
