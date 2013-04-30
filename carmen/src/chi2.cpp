#include <math.h>


double calc_chi2(double o_a, double o_b, double o_c, double o_d){
	// OBSERVED
	//    o_a  o_b
	//    o_c  o_d
	double sum_top = o_a + o_b;
	double sum_bot = o_c + o_d;
	double sum_left = o_a + o_c;
	double sum_right = o_b + o_d;
	double sum = o_a + o_b + o_c + o_d;

	double e_a = sum_top * sum_left / sum;
	double e_b = sum_top * sum_right / sum;
	double e_c = sum_bot * sum_left / sum;
	double e_d = sum_bot * sum_right / sum;
	
	double s1 = (o_a-e_a)*(o_a-e_a) / e_a;
	double s2 = (o_b-e_b)*(o_b-e_b) / e_b;
	double s3 = (o_c-e_c)*(o_c-e_c) / e_c;
	double s4 = (o_d-e_d)*(o_d-e_d) / e_d;
	return s1 + s2 + s3 + s4;
}


// below code copied from the invaluable http://www.stat.vt.edu/~sundar/java/
// another chi2 that cuts out at 10-9 is http://www.netlib.org/a/perlman
// * * Copyright (c) 2000 by Sundar Dorai-Raj
// * * @author Sundar Dorai-Raj
// * * Email: sdoraira@vt.edu
// * * This program is free software; you can redistribute it and/or
// * * modify it under the terms of the GNU General Public License 
// * * as published by the Free Software Foundation; either version 2 
// * * of the License, or (at your option) any later version, 
// * * provided that any use properly credits the author. 
// * * This program is distributed in the hope that it will be useful,
// * * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// * * GNU General Public License for more details at http://www.gnu.org 


double lngamma(double c) {
  double cof[6];
  cof[0]=76.18009172947146;
  cof[1]=-86.50532032941677;
  cof[2]=24.01409824083091;
  cof[3]=-1.231739572450155;
  cof[4]=0.01208650973866179;
  cof[5]=-0.00005395239384953e-5;
  double xx=c;
  double yy=c;
  double tmp = xx + 5.5 - (xx + 0.5) * log(xx + 5.5);
  double ser = 1.000000000190015;
  for(int j=0;j<=5;j++)
    ser += (cof[j] / ++yy);
  return(log(2.5066282746310005*ser/xx)-tmp);
}


double pchisq(double q, double df) {
  // Posten, H. (1989) American Statistician 43 p. 261-265
	double df2=df/2;
	double q2=q/2;
	double CFL, CFU, tk, p;
	int nn=5;
	if(q<=0 || df<=0) 
		return -1;
	if(q<df) {
		tk = q2 * (1-nn-df2) / (df2+2*nn-1+nn*q2 / (df2+2*nn) );
		for(int kk=nn-1;kk>1;kk--){
			tk=q2*(1-kk-df2)/(df2+2*kk-1+kk*q2/(df2+2*kk+tk));
			CFL= 1 - ( q2 / (df2+1+q2 / (df2+2+tk) ) );
			p=exp(df2*log(q2)-q2-lngamma(df2+1)-log(CFL));
		}
	}
	else{
		tk=(nn-df2)/(q2+nn);
		for(int kk=nn-1;kk>1;kk--){
			tk=(kk-df2)/(q2+kk/(1+tk));
			CFU=1+(1-df2)/(q2+1/(1+tk));
			p =1 - exp((df2-1) * log(q2)-q2-lngamma(df2)-log(CFU));
		}
	}
  return 1-p;
}


double pnorm(double z, bool upper) {
    // Algorithm AS66 Applied Statistics (1973) vol22 no.3
    // Computes P(Z<z)
    
    double ltone=7.0;
    double utzero=18.66;
    double con=1.28;
    double a1 = 0.398942280444;
    double a2 = 0.399903438504;
    double a3 = 5.75885480458;
    double a4 =29.8213557808;
    double a5 = 2.62433121679;
    double a6 =48.6959930692;
    double a7 = 5.92885724438;
    double b1 =0.398942280385;
    double b2 =3.8052e-8;
    double b3 =1.00000615302;
    double b4 =3.98064794e-4;
    double b5 =1.986153813664;
    double b6 =0.151679116635;
    double b7 =5.29330324926;
    double b8 =4.8385912808;
    double b9 =15.1508972451;
    double b10=0.742380924027;
    double b11=30.789933034;
	double b12=3.99019417011;
	double alnorm;
    if(z<0) {
		upper=!upper;
		z=-z;
    }
    if(z<=ltone || upper && z<=utzero) {
		double y = 0.5 * z * z;
		if(z>con) {
			alnorm= b1 * exp(-y) / (z-b2+b3 / (z+b4+b5 / (z-b6+b7 / (z+b8-b9 / (z+b10+b11 / (z+b12))))));
		}
		else {
			alnorm=0.5-z*(a1-a2*y/(y+a3-a4/(y+a5+a6/(y+a7))));
		}
    }
    else {
		alnorm=0;
    }
    if(!upper) 
		alnorm = 1-alnorm;
    return(alnorm);
}

double lnbeta(double a, double b) {
  return( lngamma(a) + lngamma(b) - lngamma(a+b));
}

double pbeta_raw(double x, double pin, double qin, bool lower_tail) {
  // Bosten and Battiste (1974).
  // Remark on Algorithm 179, CACM 17, p153, (1974).
	double eps=0.5*1e-8;
	double sml=1e-5;
	double lneps=log(eps);
	double lnsml=log(sml);
	bool swap_tail;
	double y=0;
	double q=0; 
	double p=0;
	double ans=0;
	double xb=0;
	double ps=0;
	double term=0;
	double n=0;
	double xi=0;
	double ib=0;
	double c=0;
	double p1=0;
	double finsum=0;
	if( (pin / (pin+qin)) <x) {
		swap_tail=true;
		y=1-x;
		p=pin;
		q=qin;
	}
	else {
		swap_tail=false;
		y=x;
		p=pin;
		q=qin;
	}
	double max;
	if( ( (p+q)*y/(p+1) ) <eps) {
		ans=0;
		
		y>sml ? max=y : max = sml;
		xb=p*log(max)-log(p)-lnbeta(p, q);
		if(xb>lnsml && y!=0)
			ans=exp(xb);
		if(swap_tail==lower_tail)
			ans=1-ans;
	}
	else {
		ps=q-floor(q);
		if(ps==0) 
			ps=1;
		xb=p*log(y)-lnbeta(ps,p)-log(p);
		ans=0;
		if(xb>=lnsml) {
			ans=exp(xb);
			term=ans*p;
			if(ps!=1) {
				lneps/log(y)>4.0 ? n = lneps/log(y) : n=4.0;
				for(int i=1;i<=n;i++) {
					xi=i;
					term*=(xi-ps)*y/xi;
					ans+=term/(p+xi);
				}
			}
		}
		if(q>1) {
			xb=p*log(y)+q*log(1-y)-lnbeta(p,q)-log(q);
			xb/lnsml>0.0 ? ib= xb/lnsml : ib=0;
			term=exp(xb-ib*lnsml);
			c=1/(1-y);
			p1=q*c/(p+q-1);
			finsum=0;
			n=q;
			if(q==n) 
				n--;
			for(int i=1;i<=n;i++) {
				if(p1<=1 && term/eps<=finsum) 
					break;
				xi=i;
				term=(q-xi+1)*c*term/(p+q-xi);
				if(term>1) {
					ib--;
					term*=sml;
				}
				if(ib==0)
					finsum += term;
			}
			ans+=finsum;
		}
		if(swap_tail==lower_tail)
			ans=1-ans;
		ans<1 ? y=ans : y = 1;
		y > 0 ? ans = y : ans=0;
	}
	return ans;
}


double pbeta(double x, double pin, double qin, bool lower_tail) {
	lower_tail=true;
	if(pin<=0 || qin<=0) 
		return -1;
	if(x<=0) 
		return (lower_tail ? 0.0:1.0);
	if(x>=1) 
		return (lower_tail ? 1.0:0.0);
	return 
		pbeta_raw(x, pin, qin, lower_tail);
}

double percentile_t(double t, double df) {
	bool lower_tail=true;
	double val;
	if( df > 400000 ) {
		val = 1.0 / (4.0 * df);
		return pnorm( t * (1.0 - val) / sqrt(1.0 + t*t*2.0*val), !lower_tail);
	}
	val = pbeta( df / (df + t*t), df / 2.0, 0.5, true);
	if(t<=0) 
		lower_tail=!lower_tail;
	val /= 2.0;
	return (lower_tail ? 1-val : val);
}	
