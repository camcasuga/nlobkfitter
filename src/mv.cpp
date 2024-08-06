#include "mv.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include "nlobk_config.hpp"
#include <iostream>
#include "solver.hpp"
#include <tools/tools.hpp>

using std::cerr; using std::endl;

double MV::DipoleAmplitude(double r, double b)
{
/*
 * cerr << "Note: MVgamma model is modified!!!" << endl;

	// rcmv from 1507.03651
	BKSolver solver;	// Becaus Alpha_s is defined there

	// GBW
	double q0sqr=0.460*0.460;
	double p=1.148;
	double n = 1.0 - std::exp(-std::pow(r*r*q0sqr/4.0,p));
	if (n<0) return 0;
	if (n>1) return 1;
	return std::pow(n, 1.0/p);
*/
	/*
	double asat = solver.Alphas(9999999)*config::NC/M_PI;
	double as = solver.Alphas(r)*config::NC/M_PI;
	const double p = 0.541;
	const double q0sqr =0.621*0.621;
	double res = 1.0 - std::exp(
		-std::pow(  r*r*q0sqr/4.0 * as *( 1.0 + std::log( asat / as ) ), p));
 	return std::pow(res, 1.0/p);
	*/
	/*
	// glauberized rcmv
	double sigma0=28.7977*2.568;
	int A=208;
	InitializeWSDistribution(A);
	double ta = T_A(anomalous_dimension, A);	// 1st parameter is called anomalous_dimension, here it is b...
	double res = 1.0 - std::exp( - std::pow(A * ta * sigma0/2.0 * r*r*qs0sqr/4.0 * as *( 1.0 + std::log( asat / as ) ), p));
///TODO NOTE: one should use different C^2 in alphas, actually C=1.550, see
//https://indico.cern.ch/event/469857/contributions/1978360/attachments/1277974/1896978/Triantafyllopoulos.pdf
// TODO~ How to exponentiate rcmv with exponents!!!!!!!!!!!!!!!!!!!!

	if (res < 0) return 0;
	if (res > 1) return 1;
	//cout << "amplitude at " << r << " is " << res << ", as " << as << " exponent " << - r*r*q0sqr/4 * as *( 1 + std::log( asat / as ) ) << endl;
	res = std::pow(res, 1.0/p);
cout << r<< " " << res << endl;
	return res;
	*/
	if (b>1e-10)
		cerr << "Impact parameter is not supported!" << LINEINFO << endl;
	if (r < 1e-30)
		return 0;
	const double e = 2.7182818;
	double exponent = std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0
		* std::log( 1.0/(r*lambdaqcd) + ec*e);
	if (exponent < 1e-5) return exponent;
	else return 1.0 - std::exp(-exponent);



void MV::SetQsqr(double qsqr)
{
	qs0sqr=qsqr;
}

void MV::SetAnomalousDimension(double gamma_)
{
	anomalous_dimension=gamma_;
}

void MV::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

void MV::SetE(double ec_)
{
	ec=ec_;
}

double MV::GetE()
{
    return ec;
}

std::string MV::GetString()
{
	std::stringstream ss;
	ss << "MV model, Q_s0^2 = " << qs0sqr << " GeV^2, \\gamma = " << anomalous_dimension
		<< ", coefficient of E inside Log is " << ec
		<< ", x0=" << x0 << ", \\Lambda_QCD = " << lambdaqcd << " GeV";
	return ss.str();
}

/*
 * Set some reasonable parameters
 * That is, MV1 model for nucleus
 */
MV::MV()
{
	qs0sqr = 0.10;
	x0=0.01;
	ec=1.0;
	lambdaqcd=config::LAMBDAQCD;
	anomalous_dimension=1;
}
