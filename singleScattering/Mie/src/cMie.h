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

#include <complex.h>

int Mie(double x, double complex m, int nt, double Theta[], double complex S1[],double complex S2[]);

int Nmax(double x, double complex m);

double complex calc_an(int n, double XL, double complex Ha, double complex mL, double complex PsiXL, double complex ZetaXL, double complex PsiXLM1, double complex ZetaXLM1);

double complex calc_bn(int n, double XL, double complex Hb, double complex mL, double complex PsiXL, double complex ZetaXL, double complex PsiXLM1, double complex ZetaXLM1);

double complex calc_S1_n(int n, double complex an, double complex bn, double Pin, double Taun);

double complex calc_S2_n(int n, double complex an, double complex bn, double Pin, double Taun);

int calc_Mie(double x, double complex m, int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double complex S1[], double complex S2[]);

void calc_scatt_prop();