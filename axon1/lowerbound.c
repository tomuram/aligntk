// now find lower-bound on distance to valid pixel in reference
if (iyv < 0)
  {
    if (ixv < 0)
      lb = hypot(hypot((double) ixv, (double) iyv),
		 IMAGE(dist, rpbw, 0, 0));
    else if (ixv >= rpbw)
      lb = hypot(hypot((double) (ixv - rpbw + 1), (double) iyv),
		 IMAGE(dist, rpbw, rpbw-1, 0));
    else
      lb = hypot((double) iyv,
		 IMAGE(dist, rpbw, ixv, 0));
  }
 else if (iyv >= rpbh)
   {
     if (ixv < 0)
       lb = hypot(hypot((double) ixv, (double) (iyv - rpbh + 1)),
		  IMAGE(dist, rpbw, 0, rpbh-1));
     else if (ixv >= rpbw)
       lb = hypot(hypot((double) (ixv - rpbw + 1),
			(double) (iyv - rpbh + 1)),
		  IMAGE(dist, rpbw, rpbw-1, rpbh-1));
     else
       lb = hypot((double) (iyv - rpbh + 1),
		  IMAGE(dist, rpbw, ixv, rpbh-1));
   }
 else
   {
     if (ixv < 0)
       lb = hypot((double) ixv,
		  IMAGE(dist, rpbw, 0, iyv));
     else if (ixv >= rpbw)
       lb = hypot((double) (ixv - rpbw + 1),
		  IMAGE(dist, rpbw, rpbw-1, iyv));
     else
       lb = IMAGE(dist, rpbw, ixv, iyv);
   }
