//I read the file where there are points of the par space I visited
//I have to compute the autocorrelation (AC) between points at progressively grater distance
//until I reach a distance such that the autocorrelation is almost zero

//HOW TO RUN THE PROGRAM: ./autocorr fileName.txt
//fileName.txt is the file downloaded from SED generator (es. S50014813-VizieR.txt)

#include <stdio.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#define f 0.01    //minimum value I set for autocorrelation
                  //I want autocorr(D) to be smaller than f for each parameter
#define d_max 50000      //maximum distance at which I compute the AC

using namespace std;

//I give N in input to skip vector size computation
double getMean(vector<double> p, int N) {
   double mean = 0;
   for (int i = 0; i < N; i++) {
      mean += p[i];
   }
   return mean / N;
}

double getVar(vector<double> p, double m, int N) {
   double Var = 0;
   for (int i = 0; i < N; i++) {
      Var += pow(p[i] - m, 2);
   }
   return Var / (N-1);
}

double getAC(vector<double> p, double m, double Var, int d, int N) {
   double ac = 0;
   for (int i = 0; i+d < N; i++) {
      ac += (p[i] - m) * (p[i+d] - m);
   }
   return ac / (N-d-1) / Var;
}

int main(int argc, char** argv) {
   string fileName = argv[1];
   fileName = "points_" + fileName;
   vector<double> *p = new vector<double> [5];   //p[0] is the vector of the steps of the random varable M, p[1]...Mr, p[2]...Rin, and so on
   string line;
   ifstream leggo(fileName);
   if (!leggo.good()) {
      cout << "File name to read is not valid!" << endl;
      exit(EXIT_FAILURE);
   }
   //I save the values of the parameters steps in the 5 vectors
   while(getline(leggo, line)) {
      stringstream ss;
      ss << line;
      double val;
      ss >> val;  //the first number in each line is the Chi
      for (int i = 0; i < 5; i++) {
         ss >> val;
         p[i].push_back(val);
      }
   }
   int N = p[0].size();    //this is the number of points in the par space for all parametrs
   double *m = new double [5];  //vector of the means of the parameters p[]
   double *Var = new double [5];    //vector of the variances of the parameters p[]
   for(int i = 0; i < 5; i++) {
      m[i] = getMean(p[i], N);
      Var[i] = getVar(p[i], m[i], N);
   }
   int d;  //distance
   for (int i = 0; i < 5; i++) {
      if(Var[i] == 0) {
         cout << "Var = 0  You froze parameter " << i << endl;
         continue;
      }
      for (d = 5; d < d_max; d=d+5) {
         double AC = getAC(p[i], m[i], Var[i], d, N);
         //cout << "d = " << d << "     " << AC << endl;
         if(AC < f) {
            cout << "d = " << d << "  for parameter " << i << endl;
            break;
         }
      }
      if(d == d_max) cout << "d over " << d_max << " for parameter " << i << endl;
   }
   cout << "Now just choose the higher distance from above (Dmax) and plot your results with:" << endl
   << "python histos.py " << argv[1] << " Dmax" << endl;
   delete[] p;
   delete[] m;
   delete[] Var;
   return 0;
}
