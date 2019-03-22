// Functions needed to generate a random number from a normal standard distribution

#ifndef normal_h
#define normal_h

#define t 5.0       //I look for step values in the interval [-t*sigma, t*sigma]
                  //at x = (+-)t*sigma f is small enough: f(x=t*sigma) / f(0) = exp(-t^2/2)

                  
//normal pdf centered in zero
double f(double x) {
   return 1/sqrt(2*pi) * exp(-pow(x,2)/2);
}

//trapezoidal integration
double integrate(double a, double b) {
   double dx = (b - a) / N;
   double sum = 0;
   double f1 = f(a);
   double f2 = f(a+dx);
   for (int i = 0; i < N; i++) {
      sum += (f1 + f2)*dx/2;
      f1 = f2;
      f2 = f(a+dx +dx*(i+1));
   }
   return sum;
}
//cumulative
double F(double x) {
   //I should integrate from -inf to x but I choose to start from x = -t*sigma where f is small enough:
   //f(x=-t*sigma) / f(0) = exp(-t^2/2)
   return integrate(-t, x);
}

double find_x(double P, double* x, double* CDF) {
   int i = N/2;   //x[N/2] = 0
   if (P >= 0.5) {
      while (P > CDF[i]) {
         i++;
      }
   }
   else {
      while (P < CDF[i]) {
         i--;
      }
   }
   return (x[1]-x[0]) * (P - CDF[i-1]) / (CDF[i] - CDF[i-1]) + x[i-1];
}

double cumulative(double *x, double *CDF) {
   double P = (double) rand() / RAND_MAX;    //random number from 0 to 1
   return find_x(P, x, CDF);
}

void save(double *x, double *CDF) {
   double dx = (double) 2*t / N;
   for (int i = 0; i < N+1; i++) {
      x[i] = -t + dx*i;
      CDF[i] = F(x[i]);
   }
   return;
}

#endif
