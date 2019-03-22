// Functions needed to compute the spectrum of the Shakura-Sunyaev model

#ifndef ss_h
#define ss_h

//cgs system
#define h 6.6E-27   //Planck constant
#define c 3.0E10    //speed of light
#define G 6.7E-8    //gravitational constant
#define k 1.4E-16   //Boltzmann consant
#define s 5.7E-5    //Stefan-Boltzmann constant (sigma)
#define n_approx 0.08      //radiative efficiency, used to estimate parameters M, Mr
#define Msun 2E33    //mass of the sun
#define yr 365*24*60*60    //seconds in a year
#define A pow(10,15.25)    //constant in A8

#define pi 3.1415926535897
#define N 500

#define constant 10E-5  //constant to keep the integrand near 1

//dimensionless flux distribution (F = Ft * P)
double Ft(double x) {
   return 3/pi * pow(x,-3) * (1-pow(x,-0.5));
}

//the luminosity is an integral, function of the frequency, I have to compute numerically
//the integrand is
double f(double x, double nu, double M, double Mr, double Rin) {
   //coefficients I will use
   double Rg = G*M/(c*c);
   double n = Rg/(2 * Rin);    //radiative efficiency
   double P = pow(n, 3) * Mr*c*c/pow(Rg,2);    //F = Ft * P
   return constant * x /(exp(h*nu/k*pow(s/P/Ft(x), 0.25)) - 1);
}

//I want to perform the integration in the logarithmic space:
//to do this I simply change variable x = e^y and get the new integrand g(y)
double g(double y, double nu, double M, double Mr, double Rin) {
   return exp(y) * f(exp(y), nu, M, Mr, Rin);
}

//trapezoidal integration
double integrate(double a, double b, double nu, double M, double Mr, double Rin) {
   double dx = (b - a) / N;
   double sum = 0;
   double g1 = g(a, nu, M, Mr, Rin);
   double g2 = g(a+dx, nu, M, Mr, Rin);
   for (int i = 0; i < N-1; i++) {
      sum += (g1 + g2)*dx/2;
      g1 = g2;
      g2 = g(a+dx + dx*(i+1), nu, M, Mr, Rin);
   }
   return sum;
}

//I consider the angle theta (line of sight - disk normal)
double NuL(double nu, double M, double Mr, double Rin, double Rout, double theta) {
   double x_out = Rout / Rin;
   //note the logarithmic bounds of the integral
   return cos(theta*pi/180) * nu * 8*pow(pi*Rin/c,2)*h*pow(nu,3) * integrate(log(1), log(x_out), nu, M, Mr, Rin) / constant;
}

#endif
