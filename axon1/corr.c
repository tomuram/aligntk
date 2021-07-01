
findcorr ()
{
  for (x = 0; x <= hw; ++x)
    lim[x] = (int) floor(sqrt((hw + 0.5) * (hw + 0.5) - x * x));

  // first compute the sums for (-1, -1)
  n = 0;
  sumA = 0.0;
  sumB = 0.0;
  sumA2 = 0.0;
  sumB2 = 0.0;
  sumAB = 0.0;
  for (y = 0; y < hw; ++y)
    {
      maxX = lim[y+1] - 1;
      for (x = 0; x <= maxX; ++x)
	{
	  if (!valid[y*w+x])
	    continue;
	  av = a[y*w+x];
	  bv = b[y*w+x];
	  startSumA += av;
	  startSumA2 += av * av;
	  startSumB += bv;
	  startSumB2 += bv * bv;
	  startSumAB += av * bv;
	  ++startN;
	}
    }
  // go through y
  for (yc = 0; yc < h; ++yc)
    {
      for (i = 1; i <= hw; ++i)
	{
	  // subtract old top
	  x = i - 1;
	  y = yc - lim[i] - 1;
	  if (y >= 0 && y < h && x < w && valid[y*w+x])
	    {
	      av = a[y*w+x];
	      bv = b[y*w+x];
	      startSumA -= av;
	      startSumA2 -= av * av;
	      startSumB -= bv;
	      startSumB2 -= bv * bv;
	      startSumAB -= av * bv;
	      --startN;
	    }

	  // add new bottom
	  y = yc + lim[i];
	  if (y < h && x < w && valid[y*w+x])
	    {
	      av = a[y*w+x];
	      bv = b[y*w+x];
	      startSumA += av;
	      startSumA2 += av * av;
	      startSumB += bv;
	      startSumB2 += bv * bv;
	      startSumAB += av * bv;
	      ++startN;
	    }
	}

      sumA = startSumA;
      sumA2 = startSumA2;
      sumB = startSumB;
      sumB2 = startSumB2;
      sumAB = startSumAB;
      n = startN;

      for (xc = 0; xc < w; ++xc)
	{
	  // subtract old left
	  x = xc - lim[0] - 1;
	  if (x >= 0 && x < w && valid[yc*w+x])
	    {
	      av = a[yc*w+x];
	      bv = b[yc*w+x];
	      sumA -= av;
	      sumA2 -= av * av;
	      sumB -= bv;
	      sumB2 -= bv * bv;
	      sumAB -= av * bv;
	      --n;
	    }

	  // add new right
	  x = xc + lim[0];
	  if (x >= 0 && x < w && valid[yc*w+x])
	    {
	      av = a[yc*w+x];
	      bv = b[yc*w+x];
	      sumA -= av;
	      sumA2 -= av * av;
	      sumB -= bv;
	      sumB2 -= bv * bv;
	      sumAB -= av * bv;
	      ++n;
	    }

	  for (i = 1; i <= hw; ++i)
	    {
	      // subtract old left
	      x = xc - lim[i] - 1;
	      if (x < 0)
		continue;
	      y = yc - i;
	      if (y >= 0 && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA -= av;
		  sumA2 -= av * av;
		  sumB -= bv;
		  sumB2 -= bv * bv;
		  sumAB -= av * bv;
		  --n;
		}
	      y = yc + i;
	      if (y < h && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA -= av;
		  sumA2 -= av * av;
		  sumB -= bv;
		  sumB2 -= bv * bv;
		  sumAB -= av * bv;
		  --n;
		}

	      // add new right
	      x = xc + lim[i];
	      if (x >= h)
		continue;
	      y = yc - i;
	      if (y >= 0 && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA += av;
		  sumA2 += av * av;
		  sumB += bv;
		  sumB2 += bv * bv;
		  sumAB += av * bv;
		  ++n;
		}
	      y = yc + i;
	      if (y < h && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA += av;
		  sumA2 += av * av;
		  sumB += bv;
		  sumB2 += bv * bv;
		  sumAB += av * bv;
		  ++n;
		}
	    }

	  // calculate the correlation for the disk surrounding (xc, yc)
	  if (n > 0)
	    {
	      meanA = sumA / n;
	      meanB = sumB / n;
	      denom = (sumA2 - 2.0 * meanA * sumA + n * meanA * meanA) *
		(sumB2 - 2.0 * meanB * sumB + n * meanB * meanB);
	      if (denom < 0.001)
		corr = -1000000.0;
	      else
		corr = (sumAB - meanA * sumB - meanB * sumA + n * meanA * meanB) / sqrt(denom);
	    }
	  else
	    corr = -1000000.0;
	}
    }
}
