float ran_pos(float);
float ran_uniform(const float, const float);
float ran_gauss(const float, const float);
float ran_gauss_any(const float, const float);
unsigned int ran_poisson (float );
  
float ran_pos(float range)
{
  float r;
  do{
    r = ran(range);
  }
  while(r == 0);
  
  return r;
}
//Unifrom with limits
float ran_uniform(const float min, const float max)
{
  return (min + ran_pos(min - max));
}

//Polar (Box-Mueller) method, from GSL, from Knuth v2, 3rd ed, p122.
float ran_gauss(const float mu, const float sigma)
{//Generate Gaussian variate with mean MU and std SIGMA, only >0 outcomes
  float x, y, r2, v = -1.0;
  while(v < 0){//Simple workaround for negative outputs
    do{
      x = -1 + 2 * ran_pos(1.0);
      y = -1 + 2 * ran_pos(1.0);
      r2 = x * x + y * y;
    }
    while(r2 > 1.0 || r2 == 0);
  
    v = mu + sigma * y * sqrt(-2.0 * log(r2) / r2);
  }
  
  return v;
}

float ran_gauss_any(const float mu, const float sigma)
{//Generate Gaussian variate with mean MU and std SIGMA
  float x, y, r2, v = -1.0;
  do{
    x = -1 + 2 * ran_pos(1.0);
    y = -1 + 2 * ran_pos(1.0);
    r2 = x * x + y * y;
  }
  while(r2 > 1.0 || r2 == 0);
  
  v = mu + sigma * y * sqrt(-2.0 * log(r2) / r2);
  
  return v;
}


// POISSON related functions. GSL CODE!
/* The poisson distribution has the form

   p(n) = (mu^n / n!) exp(-mu) 

   for n = 0, 1, 2, ... . The method used here is the one from Knuth. */
unsigned int ran_poisson (float mu)
{
  float emu;
  float prod = 1.0;
  unsigned int k = 0;

  //For large mu's (>10), not used here
  // while (mu > 10)
  //   {
  //     unsigned int m = mu * (7.0 / 8.0);

  //     double X = ran_gamma_int (m);

  //     if (X >= mu)
  //       {
  //         return k + ran_binomial (mu / X, m - 1);
  //       }
  //     else
  //       {
  //         k += m;
  //         mu -= X;
  //       }
  //   }

  /* This following method works well when mu is small */
  // For mu <= 10
  emu = exp (-mu);

  do
    {
      prod *= ran (1.0);
      k++;
    }
  while (prod > emu);

  return k - 1;

}
