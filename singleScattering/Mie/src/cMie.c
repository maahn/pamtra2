/*
    Copyright (C) 2017 - 2018 Davide Ori dori@uni-koeln.de
    Institute for Geophysics and Meteorology - University of Cologne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "cMie.h"
#include <stdio.h>
#include <math.h>

#define MAXTHETA 800

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a < _b ? _a : _b; })

void Mie(double x, double complex m) {
//    int N = Nmax(x,m);
    double Theta[MAXTHETA];
    double complex S1[MAXTHETA], S2[MAXTHETA];
    double Qext, Qabs, Qsca, Qbk;
    int nt = 0;
    int Nmax = calc_Mie(x, m, nt, Theta, &Qext, &Qsca, &Qabs, &Qbk, S1, S2);
}

int Nmax(double x, double complex m) {
// Simpler (conservative) version respect to what used by Pena(2009)
    int Nstop = round(x + 4.0*pow(x,1.0/3.0)+2);
    int xml = round(cabs(x*m));
    return max(xml,Nstop) + 15;
}

// Calculate an - equation (5)
double complex calc_an(int n, double XL, double complex Ha, double complex mL, double complex PsiXL, double complex ZetaXL, double complex PsiXLM1, double complex ZetaXLM1) {
    double complex Numer = (Ha/mL + n/XL)*PsiXL  - PsiXLM1;
    double complex Denom = (Ha/mL + n/XL)*ZetaXL - ZetaXLM1;
    return Numer/Denom;
}

// Calculate bn - equation (6)
double complex calc_bn(int n, double XL, double complex Hb, double complex mL, double complex PsiXL, double complex ZetaXL, double complex PsiXLM1, double complex ZetaXLM1) {
    double complex Numer = (Hb*mL + n/XL)*PsiXL  - PsiXLM1;
    double complex Denom = (Hb*mL + n/XL)*ZetaXL - ZetaXLM1;
    return Numer/Denom;
}

// Calculates S1_n - equation (25a)
double complex calc_S1_n(int n, double complex an, double complex bn, double Pin, double Taun) {
    return ((double)(n + n + 1)/(double)(n*n + n)) * (Pin*an + Taun*bn);
}

// Calculates S2_n - equation (25b) just switches Pin and Taun)
double complex calc_S2_n(int n, double complex an, double complex bn, double Pin, double Taun) {
    return ((double)(n + n + 1)/(double)(n*n + n)) * (Taun*an + Pin*bn);
}

int calc_Mie(double x, double complex m, int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double complex S1[], double complex S2[]) {
    int maxN = Nmax(x, m);
    double complex an, bn, anP1, bnP1, Qbktmp;
    double complex D1_lmlx[maxN + 2][2], D3_lmlx[maxN + 1][2];
    double complex D1XL[maxN + 2], D3XL[maxN + 1];
    double complex PsiZeta_lmlx[maxN + 1][2];
    double complex PsiXL[maxN + 1], ZetaXL[maxN + 1], PsiZetaXL[maxN + 1];
    double complex Ha[maxN + 1][2], Hb[maxN + 1][2];
    double Pi[maxN + 1][nTheta], Tau[maxN + 1][nTheta];
    double complex z1;
    double x2;
    int n, t;

    // Initialize the scattering parameters
    *Qext = 0;
    *Qsca = 0;
    *Qabs = 0;
    *Qbk = 0;
    Qbktmp = 0.0*I;

    // Initialize Pi, Tau and the scattering amplitudes
    for (t=0; t<nTheta; t++) {
        printf("Init Pi, Tau, S1 and S2 arrays ntheta \n");
        Pi[0][t] = 0.0;
        Tau[0][t] = 0.0;
        S1[t] = 0.0*I;
        S2[t] = 0.0*I;
    }
    
    //********************************************************//
    // Calculate D1, D3 and PsiZeta for z1 in the first layer //
    //********************************************************//
    z1 = x*m;
    // Downward recurrence for D1 - equations (16a) and (16b)
    D1_lmlx[maxN+1][1] = 0.0*I;
    for (n=maxN+1; n>0; n--) {
        D1_lmlx[n-1][1] = (double)n/z1 - 1./(D1_lmlx[n][1] + (double)n/z1);
    }

    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_lmlx[0][1] = 0.5*(1. - (cos(2.*creal(z1))+I*sin(2.*creal(z1)))*exp(-2.*cimag(z1)));
    D3_lmlx[0][1] = I;
    for (n=1; n<=maxN; n++) {
        PsiZeta_lmlx[n][1] = PsiZeta_lmlx[n-1][1]*((double)n/z1 - D1_lmlx[n-1][1])*((double)n/z1 - D3_lmlx[n-1][1]);
        D3_lmlx[n][1] = D1_lmlx[n][1] + I/PsiZeta_lmlx[n][1];
    }

  //******************************************************************//
  // Calculate Ha and Hb in the first layer - equations (7a) and (8a) //
  //******************************************************************//
    for (n=1; n<=maxN; n++) {
      Ha[n][1] = D1_lmlx[n][1];
      Hb[n][1] = D1_lmlx[n][1];
    }

  ////////////////////////////////////////
  //Calculate D1, D3 and PsiZeta for XL //
  ////////////////////////////////////////
    z1 = x + 0.0*I;

  // Downward recurrence for D1XL - equations (16a) and (16b)
    D1XL[maxN+1] = 0.0*I;
    for (n=maxN+1; n>0; n--) {
        D1XL[n-1] = (double)n/z1 - 1.0/(D1XL[n] + (double)n/z1) ;
    }

    //Upward recurrence for PsiXL, ZetaXL and D3XL - equations (18b), (18d) and (20a) - (21b)
    PsiXL[0] =  sin(creal(z1)) + I*0.0;
    ZetaXL[0] = sin(creal(z1)) - I*cos(creal(z1));

    PsiZetaXL[0] = 0.5*(1.0 - ( cos(2*creal(z1)) + I*sin(2*creal(z1)) )*exp(-2.0*cimag(z1)));
    D3XL[0] = I;
    for (n=1; n<=maxN; n++) {
        PsiXL[n]  =  PsiXL[n-1]*((double)n/z1 - D1XL[n-1]);
        ZetaXL[n] = ZetaXL[n-1]*((double)n/z1 - D3XL[n-1]);

        PsiZetaXL[n] = PsiZetaXL[n-1]*((double)n/z1 - D1XL[n-1])*((double)n/z1 - D3XL[n-1]);
        D3XL[n] = D1XL[n] + I/PsiZetaXL[n];
    }

  //*********************************************************************//
  // Finally, we calculate the an and bn coefficients and the resulting  //
  // scattering parameters                                               //
  //*********************************************************************//

    x2 = x*x;
    anP1 = calc_an(1, x, Ha[1][1], m, PsiXL[1], ZetaXL[1], PsiXL[0], ZetaXL[0]);
    bnP1 = calc_bn(1, x, Hb[1][1], m, PsiXL[1], ZetaXL[1], PsiXL[0], ZetaXL[0]);

    for (n=1; n<=maxN; n++) {
        an = anP1;
        bn = bnP1;

        anP1 = calc_an(n+1, x, Ha[n+1][1], m, PsiXL[n+1], ZetaXL[n+1], PsiXL[n], ZetaXL[n]);
        bnP1 = calc_bn(n+1, x, Hb[n+1][1], m, PsiXL[n+1], ZetaXL[n+1], PsiXL[n], ZetaXL[n]);

    // Equation (27)
        *Qext = *Qext + (double)(2*n + 1)*(creal(an) + creal(bn));
    // Equation (28)
        *Qsca = *Qsca + (double)(2*n + 1)*(creal(an)*creal(an) + cimag(an)*cimag(an) + creal(bn)*creal(bn) + cimag(bn)*cimag(bn));
    // Equation (33) inner sum
        Qbktmp = Qbktmp + (double)((2*n + 1)*(1 - 2*(n % 2))) * (an - bn);

    //****************************************************//
    // Calculate Pi_n and Tau_n for all values of Theta   //
    // Equations (26a) - (26c)                            //
    //****************************************************//
        for (t=0; t<nTheta; t++) {
            Pi[n][t] = ( (n==1) ? 1.0 : (((double)(2*n - 1)*cos(Theta[t])*Pi[n-1][t] - (double)n*Pi[n-2][t])/((double)(n-1))));
            Tau[n][t] = (double)n*cos(Theta[t])*Pi[n][t] - (double)(n+1)*Pi[n-1][t];

            S1[t] = S1[t] + calc_S1_n(n, an, bn, Pi[n][t], Tau[n][t]);
            S2[t] = S2[t] + calc_S2_n(n, an, bn, Pi[n][t], Tau[n][t]);
        }
    }

    *Qext = 2*(*Qext)/x2;                                 // Equation (27)
    *Qsca = 2*(*Qsca)/x2;                                 // Equation (28)
    *Qabs = *Qext - *Qsca;                                // Equation (30)
    *Qbk = (creal(Qbktmp)*creal(Qbktmp) + cimag(Qbktmp)*cimag(Qbktmp))/x2;    // Equation (33) last norm and xl
    printf("%+.5e, %+.5e, %+.5e, %+.5e\n", *Qext, *Qsca, *Qabs, *Qbk);
    return maxN;
}
