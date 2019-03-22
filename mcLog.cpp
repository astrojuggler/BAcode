//I read data in data_fileName.txt and do a random walk in the logarithmic parameters' space
//I write a file points_filename.txt that has 6 columns: Chi^2, log(M), log(Mr), log(Rin), log(Rout)

//HOW TO RUN THE PROGRAM: ./mcLog fileName.txt
//fileName.txt is the file downloaded from SED generator

#include <stdio.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "ss.h"
#include "normal.h"

#define pi 3.1415926535897
#define N_steps 10000    //number of steps in the random walk
#define coeff 1      //normalization of the Chi (in order to compute transitions probability)

using namespace std;

//read data and fill vectors nu, L, L_err
void getData(string fileName, vector<double>& nu, vector<double>& L, vector<double>& L_err) {
   ifstream leggo("data/"+fileName);
   if (!leggo.good()) {
      cout << "Reading file name is not valid!" << endl;
      exit(EXIT_FAILURE);
   }
   string line;
   while(getline(leggo, line)) {
      stringstream ss;
      ss << line;
      double val;
      ss >> val;
      nu.push_back(val);
      ss >> val;
      L.push_back(val);
      ss >> val;
      L_err.push_back(val);
      ss.clear();
   }
   leggo.close();
   return;
}

//read data saved in the vectors and -given a set of parametrs *p- evaluate the Chi Squared
double getChi(vector<double> nu, vector<double> L, vector<double> L_err, double *p) {
   double Chi = 0;
   double size = nu.size();
   for(int i = 0; i < size; i++) {
      // if you fix Rin, put Rin instead of pow(19,p[2])
      double m = NuL(nu[i], pow(10,p[0]), pow(10,p[1]), pow(10,p[2]), pow(10,p[3]), p[4]);
      Chi += pow((L[i] - m)/L_err[i], 2);
   }
   return Chi;
}

//propose a new point p_new making for each parameter a normal random step with sigmas sig[i]
void do_step(double *p_new, double *p, double* x, double* CDF, double *sig) {
   for (int i = 0; i < 5; i++) {
      p_new[i] = p[i] + cumulative(x, CDF) * sig[i];   //N(0,1) * s ~ N(0,s^2)
   }
   //if I want Rin to be 6*G*M/(c*c) --> p2 = log(6GM/c^2) = log(M) + log(6G/c^2)
   //p_new[2] = p_new[0] + log10(6*G/(c*c));
   //Rin is suspected to be larger: Rin = 7*2MG/c^2 look at RinFree-analysis for the hint
   //p_new[2] = p_new[0] + log10(14*G/(c*c));
   // p3 = log(Rout) = log(1000*Rin) = log(1000) + log(Rin) = 3 + p2
   //p_new[3] = 3 + p_new[2];
   return;
}

//p = p_new
void accept(double *p, double *p_new){
   for (int i = 0; i < 5; i++) {
      p[i] = p_new[i];
   }
   return;
}

//return a string that is the line to add to the file where I save the points of the random walk
string newLine(double Chi, double *p) {
   string line;
   stringstream ss;
   ss << Chi;
   for (int j = 0; j < 5; j++) {
      ss << " " << p[j];
   }
   line = ss.str();
   ss.clear();
   return line;
}

//read data in the vectors and estimate M and Mr (see estimate_par.cpp), give reasonable values to Rin Rout theta and to sigmas of all par
void setPar(vector<double> nu, vector<double> L, vector<double> L_err, double *p, double *sig) {
   double nu_p = 0, L_p = 0;    //frequency of the peak, peak luminositiy
   for (int i = 0; i < nu.size(); i++) {
      if (L[i] > L_p) {
         nu_p = nu[i];
         L_p = L[i];
      }
   }
   double Mr = 2*L_p/(n_approx*c*c);   //10E28
   double M = 1E9*Msun*pow(A/nu_p,2) * pow(n_approx/0.1,1.5) * pow(Mr*yr/Msun,0.5);     //10E44
   double Rin = 6*G*M/(c*c);  //10E17
   //double Rin = 14*G*M/(c*c);
   double Rout = 1000*Rin;    //10E20
   double theta = 30;
   p[0] = log10(M);
   //p[0] = 43.4;
   p[1] = log10(Mr);
   //p[1] = 27.8;
   p[2] = log10(Rin);
   //p[2] = 16.1;
   p[3] = log10(Rout);
   p[4] = theta;
   //I choose the sigmas to be:
   sig[0] = 0.05;
   sig[1] = 0.05;
   sig[2] = 0.03;
   sig[3] = 0.03;
   sig[4] = 0.5;
   return;
}

int main(int argc, char** argv) {
   srand(time(NULL));
   double first = rand();     //I don't know why but the first random number generated is always near 0.1 * RAND_MAX so I throw it away
   clock_t tm;
   int start = clock();
   //I save values (x, CDF) for cumulative method to generate normal random numbers
   double *x = new double [N+1];
   double *CDF = new double [N+1];
   save(x, CDF);

   string fileName = argv[1];
   string fileName1 = "data_" + fileName;
   vector<double> nu;        //frequencies of my data
   vector<double> L;         //luminosities of my data
   vector<double> L_err;     //luminosities errors of my data
   getData(fileName1, nu, L, L_err);    //I fill the vectors with my data
   double *p = new double [5];   //parameters
   double *sig = new double [5];   //sigma of the parameters

   setPar(nu, L, L_err, p, sig);  //first point of the random walk
   double Chi = getChi(nu, L, L_err, p);
   string fileName2 = "points_" + fileName;
   ofstream scrivo(fileName2);  //file where I write the points of the par space where I go: Chi, parameters
   scrivo << newLine(Chi, p) << endl;
   double *p_new = new double [5];  //possible next point
   int New = 0;
   int Same = 0;
   //random walk
   for(int i = 0; i < N_steps; i++) {
      do_step(p_new, p, x, CDF, sig);    //propose a new point
      double Chi_new = getChi(nu, L, L_err, p_new);
      if (Chi_new < Chi) {
         accept(p, p_new); //p = p_new
         Chi = Chi_new;
         scrivo << newLine(Chi_new, p_new) << endl;
         cout << "new" << endl;
         New++;
      }
      else {
         double P = (double) exp(coeff * (Chi - Chi_new));   //probability whitin I accept the new point
         cout << "coeff * (Chi - Chi_new) = " << coeff * (Chi - Chi_new) << endl;
         cout << "P = " << P << endl;
         double P_RAND = (double) rand() / RAND_MAX;
         if (P_RAND < P) {
            accept(p, p_new);  //p = p_new
            Chi = Chi_new;
            scrivo << newLine(Chi_new, p_new) << endl;
            cout << "new" << endl;
            New++;
         }
         else {
            scrivo << newLine(Chi, p) << endl;
            cout << "same" << endl;
            Same++;
         }
      }
   }
   int end = clock();
   cout << "Running time: " << (double)(end - start)/CLOCKS_PER_SEC << " seconds" << endl;
   cout << "Same/New ratio: " << (double) Same / New << endl;
   cout << fileName2 << " has been created." << endl
   << "There you find the points of the walk, let's study the autocorrelation with this command:" << endl
   << "./autocorr " << argv[1] << endl;
   scrivo.close();
   delete[] p;
   delete[] sig;
   delete[] x;
   delete[] CDF;
   return 0;
}
