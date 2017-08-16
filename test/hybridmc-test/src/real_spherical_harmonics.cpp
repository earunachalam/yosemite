/* 
   These functions are the real spherical harmonics as specified in:

   J. Mol. Struct. (THEOCHEM) 419, 19 (1997) 
   
   See also: 

   J. Chem. Phys. 136, 054501 (2012)

*/ 

#include "math.h"

const double PI=acos(-1.0);

/* These are the prefactors for the spherical harmonics. I have
   incorporated the factor (-1)^m into here. Note that m=0 is not
   multiplied by sqrt(2) */

/* Rank 2 prefactors */
const double N22sqrt2 =  0.25*sqrt(5.0/(3.0*PI)); 
const double N21sqrt2 = -0.50*sqrt(5.0/(3.0*PI));
const double N20      =  0.50*sqrt(5.0/PI);

/* Rank 3 prefactors */
const double N33sqrt2 = -sqrt(7.0/(10.0*PI))/12.0;
const double N32sqrt2 =  sqrt(7.0/(15.0*PI))/4.0;
const double N31sqrt2 = -sqrt(7.0/(6.0*PI))/2.0;
const double N30      =  sqrt(7.0/PI)/2.0;

/* Rank 4 prefactors */
const double N44sqrt2 =  sqrt(1.0/(35.0*PI))/16.0;
const double N43sqrt2 = -sqrt(1.0/(70.0*PI))/4.0;
const double N42sqrt2 =  sqrt(1.0/(5.0*PI))/4.0;
const double N41sqrt2 = -1.50*sqrt(1.0/(10.0*PI));
const double N40      =  1.50*sqrt(1.0/PI);

/* Rank 5 prefactors */
const double N55sqrt2 = -sqrt(11.0/(14*PI))/720.0;
const double N54sqrt2 =  sqrt(11.0/(35.0*PI))/144.0;
const double N53sqrt2 = -sqrt(11.0/(70.0*PI))/24.0;
const double N52sqrt2 =  sqrt(11.0/(105.0*PI))/4.0;
const double N51sqrt2 = -sqrt(11.0/(15.0*PI))/2.0;
const double N50      =  sqrt(11.0/PI)/2.0;

/* Rank 6 prefactors */
const double N66sqrt2 =  sqrt(13.0/(462.0*PI))/1440.0;
const double N65sqrt2 = -sqrt(13.0/(154.0*PI))/720.0;
const double N64sqrt2 =  sqrt(13.0/(7.0*PI))/720.0;
const double N63sqrt2 = -sqrt(13.0/(210.0*PI))/24.0;
const double N62sqrt2 =  sqrt(13.0/(210.0*PI))/4.0;
const double N61sqrt2 = -sqrt(13.0/(21.0*PI))/2.0;
const double N60      =  sqrt(13.0/PI)/2.0;

/* Rank 12 prefactors */
const double N12_12sqrt2 =  sqrt(1.0/(1352078.0*PI))/191600640.0; 
const double N12_11sqrt2 = -sqrt(1.0/(2028117.0*PI))/31933440.0;
const double N12_10sqrt2 =  sqrt(1.0/(176358.0*PI))/15966720.0;
const double N12_9sqrt2  = -sqrt(1.0/(323323.0*PI))/1451520.0;
const double N12_8sqrt2  =  sqrt(1.0/(138567.0*PI))/241920.0;
const double N12_7sqrt2  = -sqrt(1.0/(138567.0*PI))/24192.0;
const double N12_6sqrt2  =  sqrt(1.0/(4862.0*PI))/12096.0;
const double N12_5sqrt2  = -sqrt(1.0/(17017.0*PI))/576.0;
const double N12_4sqrt2  =  sqrt(1.0/(2002.0*PI))/144.0;
const double N12_3sqrt2  = -sqrt(1.0/(2002.0*PI))/12.0;
const double N12_2sqrt2  = 1.25*sqrt(1.0/(3003.0*PI));
const double N12_1sqrt2  = -2.5*sqrt(1.0/(78.0*PI));
const double N12_0       = 2.5*sqrt(1.0/PI);

/* Associated Legendre polynomials */

/* rank 2 */
double P2pos2(double);
double P2pos1(double);
double P2zer0(double);

/* rank 3 */
double P3pos3(double);
double P3pos2(double);
double P3pos1(double);
double P3zer0(double);

/* rank 4 */
double P4pos4(double);
double P4pos3(double);
double P4pos2(double);
double P4pos1(double);
double P4zer0(double);

/* rank 5 */
double P5pos5(double);
double P5pos4(double);
double P5pos3(double);
double P5pos2(double);
double P5pos1(double);
double P5zer0(double);

/* rank 6 */
double P6pos6(double);
double P6pos5(double);
double P6pos4(double);
double P6pos3(double);
double P6pos2(double);
double P6pos1(double);
double P6zer0(double);

/* rank 12 */
double P12pos12(double);
double P12pos11(double);
double P12pos10(double);
double P12pos9(double);
double P12pos8(double);
double P12pos7(double);
double P12pos6(double);
double P12pos5(double);
double P12pos4(double);
double P12pos3(double);
double P12pos2(double);
double P12pos1(double);
double P12zer0(double);

/* Rank 2 real spherical harmonics */
double S2neg2(double costheta, double phi)
{
  return N22sqrt2*P2pos2(costheta)*sin(2.0*phi);
}

double S2neg1(double costheta, double phi)
{
  return N21sqrt2*P2pos1(costheta)*sin(phi);
}

double S2zer0(double costheta, double phi)
{
  return N20*P2zer0(costheta);
}

double S2pos1(double costheta, double phi)
{
  return N21sqrt2*P2pos1(costheta)*cos(phi);
}

double S2pos2(double costheta, double phi)
{
  return N22sqrt2*P2pos2(costheta)*cos(2.0*phi);
}

/* Rank 3 real spherical harmonics */
double S3neg3(double costheta, double phi)
{
  return N33sqrt2*P3pos3(costheta)*sin(3.0*phi);
}

double S3neg2(double costheta, double phi)
{
  return N32sqrt2*P3pos2(costheta)*sin(2.0*phi);
}

double S3neg1(double costheta, double phi)
{
  return N31sqrt2*P3pos1(costheta)*sin(phi);
}

double S3zer0(double costheta, double phi)
{
  return N30*P3zer0(costheta);
}

double S3pos1(double costheta, double phi)
{
  return N31sqrt2*P3pos1(costheta)*cos(phi);
}

double S3pos2(double costheta, double phi)
{
  return N32sqrt2*P3pos2(costheta)*cos(2.0*phi);
}

double S3pos3(double costheta, double phi)
{
  return N33sqrt2*P3pos3(costheta)*cos(3.0*phi);
}

/* Rank 4 real spherical harmonics */
double S4neg4(double costheta, double phi)
{
  return N44sqrt2*P4pos4(costheta)*sin(4.0*phi);
}

double S4neg3(double costheta, double phi)
{
  return N43sqrt2*P4pos3(costheta)*sin(3.0*phi);
}

double S4neg2(double costheta, double phi)
{
  return N42sqrt2*P4pos2(costheta)*sin(2.0*phi);
}

double S4neg1(double costheta, double phi)
{
  return N41sqrt2*P4pos1(costheta)*sin(phi);
}

double S4zer0(double costheta, double phi)
{
  return N40*P4zer0(costheta);
}

double S4pos1(double costheta, double phi)
{
  return N41sqrt2*P4pos1(costheta)*cos(phi);
}

double S4pos2(double costheta, double phi)
{
  return N42sqrt2*P4pos2(costheta)*cos(2.0*phi);
}

double S4pos3(double costheta, double phi)
{
  return N43sqrt2*P4pos3(costheta)*cos(3.0*phi);
}

double S4pos4(double costheta, double phi)
{
  return N44sqrt2*P4pos4(costheta)*cos(4.0*phi);
}

/* Rank 5 real spherical harmonics */
double S5neg5(double costheta, double phi)
{
  return N55sqrt2*P5pos5(costheta)*sin(5.0*phi);
}

double S5neg4(double costheta, double phi)
{
  return N54sqrt2*P5pos4(costheta)*sin(4.0*phi);
}

double S5neg3(double costheta, double phi)
{
  return N53sqrt2*P5pos3(costheta)*sin(3.0*phi);
}

double S5neg2(double costheta, double phi)
{
  return N52sqrt2*P5pos2(costheta)*sin(2.0*phi);
}

double S5neg1(double costheta, double phi)
{
  return N51sqrt2*P5pos1(costheta)*sin(phi);
}

double S5zer0(double costheta, double phi)
{
  return N50*P5zer0(costheta);
}

double S5pos1(double costheta, double phi)
{
  return N51sqrt2*P5pos1(costheta)*cos(phi);
}

double S5pos2(double costheta, double phi)
{
  return N52sqrt2*P5pos2(costheta)*cos(2.0*phi);
}

double S5pos3(double costheta, double phi)
{
  return N53sqrt2*P5pos3(costheta)*cos(3.0*phi);
}

double S5pos4(double costheta, double phi)
{
  return N54sqrt2*P5pos4(costheta)*cos(4.0*phi);
}

double S5pos5(double costheta, double phi)
{
  return N55sqrt2*P5pos5(costheta)*cos(5.0*phi);
}

/* Rank 6 real spherical harmonics */
double S6neg6(double costheta, double phi)
{
  return N66sqrt2*P6pos6(costheta)*sin(6.0*phi);
}

double S6neg5(double costheta, double phi)
{
  return N65sqrt2*P6pos5(costheta)*sin(5.0*phi);
}

double S6neg4(double costheta, double phi)
{
  return N64sqrt2*P6pos4(costheta)*sin(4.0*phi);
}

double S6neg3(double costheta, double phi)
{
  return N63sqrt2*P6pos3(costheta)*sin(3.0*phi);
}

double S6neg2(double costheta, double phi)
{
  return N62sqrt2*P6pos2(costheta)*sin(2.0*phi);
}

double S6neg1(double costheta, double phi)
{
  return N61sqrt2*P6pos1(costheta)*sin(phi);
}

double S6zer0(double costheta, double phi)
{
  return N60*P6zer0(costheta);
}

double S6pos1(double costheta, double phi)
{
  return N61sqrt2*P6pos1(costheta)*cos(phi);
}

double S6pos2(double costheta, double phi)
{
  return N62sqrt2*P6pos2(costheta)*cos(2.0*phi);
}

double S6pos3(double costheta, double phi)
{
  return N63sqrt2*P6pos3(costheta)*cos(3.0*phi);
}

double S6pos4(double costheta, double phi)
{
  return N64sqrt2*P6pos4(costheta)*cos(4.0*phi);
}

double S6pos5(double costheta, double phi)
{
  return N65sqrt2*P6pos5(costheta)*cos(5.0*phi);
}

double S6pos6(double costheta, double phi)
{
  return N66sqrt2*P6pos6(costheta)*cos(6.0*phi);
}

/* Rank 12 real spherical harmonics */
double S12neg12(double costheta, double phi)
{
  return N12_12sqrt2*P12pos12(costheta)*sin(12.0*phi);
}

double S12neg11(double costheta, double phi)
{
  return N12_11sqrt2*P12pos11(costheta)*sin(11.0*phi);
}

double S12neg10(double costheta, double phi)
{
  return N12_10sqrt2*P12pos10(costheta)*sin(10.0*phi);
}

double S12neg9(double costheta, double phi)
{
  return N12_9sqrt2*P12pos9(costheta)*sin(9.0*phi);
}

double S12neg8(double costheta, double phi)
{
  return N12_8sqrt2*P12pos8(costheta)*sin(8.0*phi);
}

double S12neg7(double costheta, double phi)
{
  return N12_7sqrt2*P12pos7(costheta)*sin(7.0*phi);
}

double S12neg6(double costheta, double phi)
{
  return N12_6sqrt2*P12pos6(costheta)*sin(6.0*phi);
}

double S12neg5(double costheta, double phi)
{
  return N12_5sqrt2*P12pos5(costheta)*sin(5.0*phi);
}

double S12neg4(double costheta, double phi)
{
  return N12_4sqrt2*P12pos4(costheta)*sin(4.0*phi);
}

double S12neg3(double costheta, double phi)
{
  return N12_3sqrt2*P12pos3(costheta)*sin(3.0*phi);
}

double S12neg2(double costheta, double phi)
{
  return N12_2sqrt2*P12pos2(costheta)*sin(2.0*phi);
}

double S12neg1(double costheta, double phi)
{
  return N12_1sqrt2*P12pos1(costheta)*sin(phi);
}

double S12zer0(double costheta, double phi)
{
  return N12_0*P12zer0(costheta);
}

double S12pos1(double costheta, double phi)
{
  return N12_1sqrt2*P12pos1(costheta)*cos(phi);
}

double S12pos2(double costheta, double phi)
{
  return N12_2sqrt2*P12pos2(costheta)*cos(2.0*phi);
}

double S12pos3(double costheta, double phi)
{
  return N12_3sqrt2*P12pos3(costheta)*cos(3.0*phi);
}

double S12pos4(double costheta, double phi)
{
  return N12_4sqrt2*P12pos4(costheta)*cos(4.0*phi);
}

double S12pos5(double costheta, double phi)
{
  return N12_5sqrt2*P12pos5(costheta)*cos(5.0*phi);
}

double S12pos6(double costheta, double phi)
{
  return N12_6sqrt2*P12pos6(costheta)*cos(6.0*phi);
}

double S12pos7(double costheta, double phi)
{
  return N12_7sqrt2*P12pos7(costheta)*cos(7.0*phi);
}

double S12pos8(double costheta, double phi)
{
  return N12_8sqrt2*P12pos8(costheta)*cos(8.0*phi);
}

double S12pos9(double costheta, double phi)
{
  return N12_9sqrt2*P12pos9(costheta)*cos(9.0*phi);
}

double S12pos10(double costheta, double phi)
{
  return N12_10sqrt2*P12pos10(costheta)*cos(10.0*phi);
}

double S12pos11(double costheta, double phi)
{
  return N12_11sqrt2*P12pos11(costheta)*cos(11.0*phi);
}

double S12pos12(double costheta, double phi)
{
  return N12_12sqrt2*P12pos12(costheta)*cos(12.0*phi);
}


/* Associated Legendre polynomials. For the real spherical harmonics,
   only m >= 0 are needed. Expressions taken from Mathematica, and
   therefore include the (-1)^m Condon-Shortley phase factor in the
   definition. */

/* Rank 2 */

double P2pos2(double x)
{
  return -3.0*(-1.0 + x*x);
}

double P2pos1(double x)
{
  return -3.0*x*sqrt(1.0 - x*x);
}

double P2zer0(double x)
{
  return 0.5*(-1.0 + 3.0*x*x);
}

/* Rank 3 */
double P3pos3(double x)
{
  return -15.0*sqrt( (1.0 - x*x)*(1.0 - x*x)*(1.0 - x*x) );
}

double P3pos2(double x)
{
  return -15.0*x*(-1.0 + x*x);
}

double P3pos1(double x)
{
  return -1.5*sqrt(1.0 - x*x)*(-1.0 + 5.0*x*x);
}

double P3zer0(double x)
{
  return 0.5*(-3.0*x + 5*x*x*x);
}

/* Rank 4 */
double P4pos4(double x)
{
  return 105.0*(-1.0 + x*x)*(-1.0 + x*x);
}

double P4pos3(double x)
{
  return -105.0*x*sqrt((1.0 - x*x)*(1.0 - x*x)*(1.0 - x*x));
}

double P4pos2(double x)
{
  return -7.5*(-1.0 + x*x)*(-1.0 + 7.0*x*x);
}

double P4pos1(double x)
{
  return -2.50*sqrt(1.0 - x*x)*(-3.0*x + 7.0*x*x*x);
}

double P4zer0(double x)
{
  return 0.125*(3.0 - 30.0*x*x + 35.0*x*x*x*x);
}

/* Rank 5 */
double P5pos5(double x)
{
  return -945.0*sqrt( (1.0-x*x)*(1.0-x*x)*(1.0-x*x)*(1.0-x*x)*(1.0-x*x) );
}

double P5pos4(double x)
{
  return 945.0*x*(-1.0+x*x)*(-1.0+x*x);
}

double P5pos3(double x)
{
  return 52.5*sqrt(1.0-x*x)*(-1.0+x*x)*(-1.0+9*x*x);
}

double P5pos2(double x)
{
  return -52.5*(-1.0+x*x)*(-x+3*x*x*x);
}

double P5pos1(double x)
{
  return -1.875*sqrt(1.0-x*x)*(1.0 - 14.0*x*x + 21.0*x*x*x*x);
}

double P5zer0(double x)
{
  return 0.125*(15.0*x -70.0*x*x*x + 63.0*x*x*x*x*x);
}

/* Rank 6 */
double P6pos6(double x)
{
  return -10395.0*(-1.0+x*x)*(-1.0+x*x)*(-1.0+x*x);
}

double P6pos5(double x)
{
  return -10395.0*x*sqrt( (1.0-x*x)*(1.0-x*x)*(1.0-x*x)*(1.0-x*x)*(1.0-x*x) );
}

double P6pos4(double x)
{
  return 472.5*(-1.0+x*x)*(-1.0+x*x)*(-1.0 + 11.0*x*x);
}

double P6pos3(double x)
{
  return 157.5*sqrt(1.0 - x*x)*(-1.0 +x*x)*(-3.0*x + 11.0*x*x*x);
}

double P6pos2(double x)
{
  return -13.125*(-1.0 + x*x)*(1.0 - 18.0*x*x +33.0*x*x*x*x);
}

double P6pos1(double x)
{
  return -2.625*sqrt(1.0-x*x)*(5.0*x - 30.0*x*x*x + 33.0*x*x*x*x*x);
}

double P6zer0(double x)
{
  return 0.0625*(-5.0 + 105.0*x*x - 315.0*x*x*x*x + 231.0*x*x*x*x*x*x);
}

/* Rank 12 */
double P12pos12(double x)
{
  return 316234143225.0*pow((-1.0 + x*x),6);
}

double P12pos11(double x)
{
  return -316234143225.0*x*sqrt( pow((1.0 - x*x),11) );
}

double P12pos10(double x)
{
  return -6874655287.5*pow((-1.0 + x*x),5)*(-1.0 + 23.0*x*x);
}

double P12pos9(double x)
{
  return -2291551762.5*sqrt(1.0-x*x)*pow((-1.0 + x*x),4)*(-3.0*x + 23.0*x*x*x);
}

double P12pos8(double x)
{
  return 81841134.375*pow((-1.0 + x*x),4)*(1.0 - 42.0*x*x + 161.0*x*x*x*x);
}

double P12pos7(double x)
{
  return 16368226.875*sqrt(1.0-x*x)*pow((-1.0 + x*x),3)*(5.0*x - 70.0*x*x*x + 161.0*x*x*x*x*x);
}

double P12pos6(double x)
{
  return -143580.9375*pow((-1.0 + x*x),3)*(-5.0 + 285*x*x - 1995.0*x*x*x*x + 3059.0*x*x*x*x*x*x);
}

double P12pos5(double x)
{
  return -143580.9375*sqrt(1.0-x*x)*pow((-1.0 + x*x),2)*(-5.0*x + 95.0*x*x*x - 399.0*x*x*x*x*x + 437.0*x*x*x*x*x*x*x);
}

double P12pos4(double x)
{
  return 1055.7421875*pow((-1.0 + x*x),2)*(5.0 - 340.0*x*x + 3230*x*x*x*x - 9044.0*x*x*x*x*x*x + 7429.0*x*x*x*x*x*x*x*x);
}

double P12pos3(double x)
{
  return 117.3046875*sqrt(1.0-x*x)*(-1.0 + x*x)*(45.0*x - 1020.0*pow(x,3) + 5814.0*pow(x,5) - 11628.0*pow(x,7) + 7429.0*pow(x,9));
}

double P12pos2(double x)
{
  return -11.73046875*(-1.0 + x*x)*(-3.0 + 225.0*x*x - 2550.0*pow(x,4) + 9690.0*pow(x,6) - 14535.0*pow(x,8) + 7429.0*pow(x,10));
}

double P12pos1(double x)
{
  return -0.15234375*sqrt(1.0-x*x)*(-231.0*x + 5775.0*pow(x,3) - 39270.0*pow(x,5) + 106590.0*pow(x,7) - 124355.0*pow(x,9) + 52003.0*pow(x,11));
}

double P12zer0(double x)
{
  return (231.0 - 18018.0*x*x + 225225.0*pow(x,4) - 1021020.0*pow(x,6) + 2078505.0*pow(x,8) - 1939938.0*pow(x,10) + 676039.0*pow(x,12))/1024.0;
}
