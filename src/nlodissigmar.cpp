#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
//#include <chrono>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include <amplitudelib/amplitudelib.hpp>

#include "nlodis_config.hpp"
#include "nlodissigmar.hpp"
#include "nlodissigmar_massiveq.hpp"
#include "helper.hpp"


namespace sigmar_config{
    // double rmax=30.0; //300.0;
    // double rmin=1e-5;
    double Nc=3.0;
    double sumef=6.0/9.0; // light quarks uds only.
    double CF=4.0/3.0; // (Nc()*Nc()-1.0)/(2.0*Nc());
    double alphaem=1.0/137.0;
    double lambdaqcd=0.241; //GeV
    double structurefunfac=1./(Sq(2*M_PI)*alphaem);
}
using namespace sigmar_config;

///===========================================================================================
//// MAIN FUNCTIONS HERE

// at low Q² Z exchange negligible ( we are only modeling gamma exchange) and the reduced cross section simplifies into
// sigma_r = F_2 - f(y) F_L = F_2 - y^2 / Y_+ * F_L = F_2 - y^2/(1+(1-y)^2) * F_L
double ComputeSigmaR::SigmarLO (double Q , double xbj, double y) {
    double FL = Structf_LLO(Q,xbj);
    double FT = Structf_TLO(Q,xbj);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigmar = F2 - fy*FL;
    return sigmar;
}
double ComputeSigmaR::SigmarLOmass (double Q , double xbj, double y, bool charm) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOpMass(Q,xbj,charm);
    double FT = fac*TLOpMass(Q,xbj,charm);
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

// UNSUB ---------------------
double ComputeSigmaR::SigmarNLOunsub ( double Q , double xbj, double y) {
    double FLic  = Structf_LLO(Q,icX0);
    double FTic  = Structf_TLO(Q,icX0);
    double FLdip = Structf_LNLOdip(Q,xbj);
    double FTdip = Structf_TNLOdip(Q,xbj);
    double FLqg  = Structf_LNLOqg_unsub(Q,xbj);
    double FTqg  = Structf_TNLOqg_unsub(Q,xbj);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    // Printing intermediate F2 FL FT
    // cout << "SigmarNLOunsub(Q xbj y) sigmar F2 FL FT "  << Q << " "
    //                                                     << xbj << " "
    //                                                     << y << " "
    //                                                     << sigma << " "
    //                                                     << F2 << " "
    //                                                     << FL << " "
    //                                                     << FT << endl;
    return sigma;
}

double ComputeSigmaR::SigmarNLOunsub_massive ( double Q , double xbj, double y, double q_mass) {
    double FLic  = Structf_LLO_massive(Q, icX0, q_mass);
    double FTic  = Structf_TLO_massive(Q, icX0, q_mass);
    double FLdip = Structf_LNLOdip_massive(Q, xbj, q_mass);
    double FTdip = Structf_TNLOdip_massive(Q, xbj, q_mass);
    double FLqg  = Structf_LNLOqg_unsub_massive(Q, xbj, q_mass);
    double FTqg  = Structf_TNLOqg_unsub_massive(Q, xbj, q_mass);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;

    // Printing intermediate F2 FL FT
    // cout << "SigmarNLOunsub_massive(Q xbj y mf) sigmar F2 FL FT "  << Q << " "
    //                                                                 << xbj << " "
    //                                                                 << y << " "
    //                                                                 << q_mass << " "
    //                                                                 << sigma << " "
    //                                                                 << F2 << " "
    //                                                                 << FL << " "
    //                                                                 << FT << endl;
    return sigma;
}

double ComputeSigmaR::SigmarNLOunsub_UniformZ2Bound ( double Q , double xbj, double y) {
    double FLic  = Structf_LLO(Q,icX0);
    double FTic  = Structf_TLO(Q,icX0);
    double FLdip = Structf_LNLOdip_z2(Q,xbj);
    double FTdip = Structf_TNLOdip_z2(Q,xbj);
    double FLqg  = Structf_LNLOqg_unsub(Q,xbj);
    double FTqg  = Structf_TNLOqg_unsub(Q,xbj);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOunsub_sigma3 ( double Q , double xbj, double y) {
    double FL3  = Structf_LNLOsigma3(Q,xbj);
    double FT3  = Structf_TNLOsigma3(Q,xbj);
    double F2 = FL3+FT3;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL3;
    return sigma;
}


// SUB ---------------------
double ComputeSigmaR::SigmarNLOsub ( double Q , double xbj, double y) {
    double FLlo  = Structf_LLO(Q,xbj);
    double FTlo  = Structf_TLO(Q,xbj);
    double FLdip = Structf_LNLOdip(Q,xbj);
    double FTdip = Structf_TNLOdip(Q,xbj);
    double FLqg  = Structf_LNLOqg_sub(Q,xbj);
    double FTqg  = Structf_TNLOqg_sub(Q,xbj);
    double FL = FLlo + FLdip + FLqg;
    double FT = FTlo + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOsub_UniformZ2Bound ( double Q , double xbj, double y) {
    double FLlo  = Structf_LLO(Q,xbj);
    double FTlo  = Structf_TLO(Q,xbj);
    double FLdip = Structf_LNLOdip_z2(Q,xbj);
    double FTdip = Structf_TNLOdip_z2(Q,xbj);
    double FLqg  = Structf_LNLOqg_sub(Q,xbj);
    double FTqg  = Structf_TNLOqg_sub(Q,xbj);
    double FL = FLlo + FLdip + FLqg;
    double FT = FTlo + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

// RISTO formulations ---------------------
double ComputeSigmaR::SigmarNLOunsubRisto ( double Q , double xbj, double y) {
    double fac = structurefunfac*Sq(Q);
    double FLic  = fac*LLOp(Q,icX0);
    double FTic  = fac*TLOp(Q,icX0);
    double FLdip = fac*LNLOdip(Q,xbj);
    double FTdip = fac*TNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgunsubRisto(Q,xbj);
    double FTqg = fac*TNLOqgunsubRisto(Q,xbj);
    double FL = FLic + FLdip + FLqg;
    double FT = FTic + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

double ComputeSigmaR::SigmarNLOsubRisto ( double Q , double xbj, double y) {
    double fac = structurefunfac*Sq(Q);
    double FLlo  = fac*LLOp(Q,xbj);
    double FTlo  = fac*TLOp(Q,xbj);
    double FLdip = fac*LNLOdip(Q,xbj);
    double FTdip = fac*TNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgsubRisto(Q,xbj);
    double FTqg = fac*TNLOqgsubRisto(Q,xbj);
    double FL = FLlo + FLdip + FLqg;
    double FT = FTlo + FTdip + FTqg;
    double F2 = FL+FT;
    double fy = Sq(y)/(1+Sq(1-y));
    double sigma = F2 - fy*FL;
    return sigma;
}

// STRUCTURE FUNCTIONS -------------------------
double ComputeSigmaR::Structf_LLO ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOp( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TLO ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TLOp( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LLO_massive ( double Q , double xbj, double m_f ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LLOp_massive(Q, xbj, m_f);
    return FL;
}

double ComputeSigmaR::Structf_TLO_massive ( double Q , double xbj, double m_f ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TLOp_massive(Q, xbj, m_f);
    return FT;
}

double ComputeSigmaR::Structf_LNLOdip ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOdip( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOdip ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOdip( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LNLOdip_massive ( double Q , double xbj, double m_f ) {
    double fac = structurefunfac*Sq(Q);
    double lilogconst = LNLOdip_massive_LiLogConst(Q, xbj, m_f);
    double iab = LNLOdip_massive_Iab(Q, xbj, m_f);
    double icd = LNLOdip_massive_Icd(Q, xbj, m_f);

    // cout << "LNLOdip raw lilog iab icd: " << lilogconst << " " << iab << " " << icd << endl;
    double FL = fac*(lilogconst + iab + icd);
    return FL;
}

double ComputeSigmaR::Structf_TNLOdip_massive ( double Q , double xbj, double m_f ) {
    double fac = structurefunfac*Sq(Q);
    double ft_dip_0 = TNLOdip_massive_I1(Q, xbj, m_f);
    double ft_dip_1 = TNLOdip_massive_I2(Q, xbj, m_f);
    double ft_dip_2 = TNLOdip_massive_I3(Q, xbj, m_f);
    
    double FT = fac*(ft_dip_0 + ft_dip_1 + ft_dip_2);
    return FT;
}

double ComputeSigmaR::Structf_LNLOdip_z2 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOdip_z2( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOdip_z2 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOdip_z2( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_LNLOqg_sub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOqgsub( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_LNLOqg_unsub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOqgunsub( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_LNLOqg_unsub_massive ( double Q , double xbj, double m_f ) {
    double fac = structurefunfac*Sq(Q);
    // double FL = fac*LNLOqgunsub_massive( Q , xbj, m_f ); // compute qg term in one go
    double FL_I1 = LNLOqgunsub_massive_I1( Q , xbj, m_f );
    double FL_I2 = LNLOqgunsub_massive_I2( Q , xbj, m_f );
    double FL_I3 = LNLOqgunsub_massive_I3( Q , xbj, m_f );

    // cout << "LNLOqg raw I1 I2 I3: " << FL_I1 << " " << FL_I2 << " " << FL_I3 << endl;
    return fac*(FL_I1 + FL_I2 + FL_I3); // FL_tot
}

double ComputeSigmaR::Structf_TNLOqg_sub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOqgsub( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_TNLOqg_unsub ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOqgunsub( Q , xbj );
    return FT;
}

double ComputeSigmaR::Structf_TNLOqg_unsub_massive ( double Q , double xbj, double m_f ) {
    double fac = structurefunfac*Sq(Q);
    double sig_I1 = TNLOqgunsub_massive_I1(Q, xbj, m_f);
    double sig_I2 = TNLOqgunsub_massive_I2(Q, xbj, m_f);
    double sig_I3 = TNLOqgunsub_massive_I3(Q, xbj, m_f);
    double FT = fac*(sig_I1 + sig_I2 + sig_I3);
    return FT;
}

double ComputeSigmaR::Structf_LNLOsigma3 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FL = fac*LNLOsigma3( Q , xbj );
    return FL;
}

double ComputeSigmaR::Structf_TNLOsigma3 ( double Q , double xbj ) {
    double fac = structurefunfac*Sq(Q);
    double FT = fac*TNLOsigma3( Q , xbj );
    return FT;
}


// obsolete?
double ComputeSigmaR::Structf_LFULLNLOunsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FLic  = fac*LLOp(Q,icX0);
    //double FLdip = fac*LNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgunsub(Q,xbj);
    double FL = FLic + FLqg;
    return FL;
}

double ComputeSigmaR::Structf_TFULLNLOunsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FTic  = fac*TLOp(Q,icX0);
    //double FTdip = fac*TNLOdip(Q,xbj);
    double FTqg = fac*TNLOqgunsub(Q,xbj);
    double FT = FTic + FTqg;
    return FT;
}

double ComputeSigmaR::Structf_LFULLNLOsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FLlo  = fac*LLOp(Q,xbj);
    //double FLdip = fac*LNLOdip(Q,xbj);
    double FLqg = fac*LNLOqgsub(Q,xbj);
    double FL = FLlo + FLqg;
    return FL;
}

double ComputeSigmaR::Structf_TFULLNLOsub ( double Q , double xbj) {
    double fac = structurefunfac*Sq(Q);
    double FTlo  = fac*TLOp(Q,xbj);
    //double FTdip = fac*TNLOdip(Q,xbj);
    double FTqg = fac*TNLOqgsub(Q,xbj);
    double FT = FTlo + FTqg;
    return FT;
}


///===========================================================================================
double NLODISFitter::operator()(const std::vector<double>& par) const
{
    double chisqr = 0;
    // params to IC
    double qs0sqr       = par[ parameters.Index("qs0sqr")];
    double e_c          = par[ parameters.Index("e_c")];
    //double fitsigma0    = 2.568*par[ parameters.Index("fitsigma0")];  // 1mb = 2.568 GeV² -- unit change into GeV
    double alphas_scaling       = par[ parameters.Index("alphascalingC2")]; // MATCH THIS IN IMPACTFACTOR ALPHA_S WHEN NLO
    double anomalous_dimension  = par[ parameters.Index("anomalous_dimension")];
    double icx0_nlo_impfac  = par[ parameters.Index("icx0_nlo_impfac")];
    double icx0_bk  = par[ parameters.Index("icx0_bk")];
    double initialconditionY0  = par[ parameters.Index("initialconditionY0")];
    double icTypicalPartonVirtualityQ0sqr  = par[ parameters.Index("icTypicalPartonVirtualityQ0sqr")];
    double old_sigma02 = par[ parameters.Index("old_sigma02")];
    double mass_charm = par[ parameters.Index("mass_charm")];
    double mass_bottom = par[ parameters.Index("mass_bottom")];
    double qMass_light  = 0.14; // GeV --- doesn't improve fit at LO
    double qMass_u = 0.00216; // GeV, PDG MSbar value
    double qMass_d = 0.00467; // GeV, PDG MSbar value
    double qMass_s = 0.093; // GeV, PDG MSbar value
    // double qMass_charm = 1.27 // PDG MSbar value // 1.35; == old default
    double qMass_charm = mass_charm;
    double qMass_b = 4.750; // GeV, pole mass scheme val
    double qMass_b_var = mass_bottom;
    bool useMasses = nlodis_config::USE_MASSES;
    bool useCharm = false;


    if (std::abs(initialconditionY0) > 1e-4)
    {
        cerr << "WARNING! initialcconditionY0" << initialconditionY0 << ", but values != 0 are not guaranteed to work! Talk to Heikki before using this!" << endl;
    }
    if (qs0sqr < 0.0001 or qs0sqr > 100 or alphas_scaling < 0.01 or alphas_scaling > 99999
        /*or fitsigma0 < 0.1 or fitsigma0 > 999 */
        or e_c < 1)
    return 9999999;

    cout << "=== Initializing Chi^2 regression === "<< " parameters (" << PrintVector(par) << ")" << endl;
    // Manual limiting for parameter range as a fail safe for Minuit2 craziness:
    /*
    if (  (qs0sqr < 0)    || (alphas_scaling < 0)     || (anomalous_dimension < 0) || (e_c < 0.4) || (initialconditionX0 < 0.01) ||
            (qs0sqr > 100)  || (alphas_scaling > 1000)  || (anomalous_dimension > 2) ) {
            chisqr = 1e7 * ( 1.0 + abs(qs0sqr) + abs(alphas_scaling) + abs(anomalous_dimension) + abs(e_c - 0.4)); // Should be some crazy large chisqr to deter Minuit2 from using this kind of parametrization.
            cout << endl << "# Calculated ILL PARAMETER chi^2 = " << chisqr  << ", parameters (" << PrintVector(par) << ")" << endl<<endl;
            return chisqr;
        }
    */

    /*
    // ***Solve resummed BK***
    */
    // cout << "=== Initialize BK solver ===" << endl;
   
    MV ic;                                            // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);                                     // e_c of MVe parametrization

    Dipole dipole(&ic);
    dipole.SetX0(icx0_bk);
    BKSolver solver(&dipole);
    double maxy = std::log(icx0_bk/(1e-5)) + initialconditionY0; // divisor=smallest HERA xbj in Q^2 range (1E-05)?
    
    double maxq2=50;     // Todo, take from actual data!
    if (useImprovedZ2Bound)
        maxy += std::log(maxq2 / icTypicalPartonVirtualityQ0sqr);

    // cout << "=== Solving BK ===" << endl;

    solver.SetAlphasScaling(alphas_scaling);
    // solver.SetEta0(par[ parameters.Index("eta0")]);
    solver.SetX0(icx0_bk);
    solver.SetICX0_nlo_impfac(icx0_nlo_impfac);
    solver.SetICTypicalPartonVirtualityQ0sqr(icTypicalPartonVirtualityQ0sqr);
    //solver.SetTmpOutput("tmp_datafile.dat");
    solver.Solve(maxy);

    // solver.GetDipole()->Save("output_dipole_uksi_balsd_heraII.dat");
    // solver.GetDipole()->Save("output_dipole_urpbsi_hera_resumpd_ic_nlobk.dat");
    //solver.GetDipole()->Save("output_dipole_lobk_x0=10_euler_KINCOST_DESTEP0.08.dat");
    //solver.GetDipole()->Save("output_dipole_lobk_x0=10_rungekutta_DESTEP0.4.dat");
    // AmplitudeLib* DipoleAmplitude_ptr; // Forward declaration of the dipole object to be initialized from a file or solved data.
    // DipoleAmplitude_ptr = new AmplitudeLib("output_dipole_uksi_balsd_heraII.dat");      // read data from existing file.
    // AmplitudeLib DipoleAmplitude(*DipoleAmplitude_ptr);

    
    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());
    DipoleAmplitude.SetInterpolationMethod(LINEAR_LINEAR);
    DipoleAmplitude.SetX0(icx0_bk);
    DipoleAmplitude.SetOutOfRangeErrors(false);
    AmplitudeLib *DipolePointer = &DipoleAmplitude;

    ComputeSigmaR SigmaComputer(DipolePointer);
    SigmaComputer.SetX0(icx0_nlo_impfac);
    SigmaComputer.SetX0_BK(icx0_bk);
    SigmaComputer.SetY0(initialconditionY0);
    SigmaComputer.SetQ0Sqr(icTypicalPartonVirtualityQ0sqr);
    SigmaComputer.SetQuarkMassLight(qMass_light);
    SigmaComputer.SetQuarkMassCharm(qMass_charm);
    SigmaComputer.SetAlphasScalingC2(alphas_scaling); // running coupling alpha_s scaling parameter C^2

    // Set running coupling and rapidity function pointters based on the setting enumerators in nlodis_config.hpp
    SigmaComputer.MetaPrescriptionSetter();

    // CUBA Monte Carlo integration library algorithm setter
    SigmaComputer.SetCubaMethod(cubaMethod);

    cout << "=== Computing Reduced Cross sections ===" << endl;

    /*
     * Loop over datapoints and compute theoretical predictions
     */
    int points=0, totalpoints = 0;
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
        totalpoints += datasets[dataset]->NumOfPoints();
    
    double fitsigma0 = 1;
    std::vector<double> datavals(totalpoints);
    std::vector<double> dataerrs(totalpoints);
    std::vector<double> thdata(totalpoints);  
    std::vector<double> var_xbj(totalpoints);
    std::vector<double> var_y(totalpoints);
    std::vector<double> var_qsqr(totalpoints);
    
    ///TODO: Does not yet fully support weight factors ,
    // after I have written a separate code that automatically finds optimal sigma02

    // These loops are trivially parallerizable
    // We only parallerize the inner loop where we have about
    // 250 points (total sigmar) and 50 points (charm)
    for (unsigned int dataset=0; dataset<datasets.size(); dataset++)
    {
#ifdef PARALLEL_CHISQR
        //reduction(+:chisqr)
    #pragma omp parallel for schedule(dynamic) reduction(+:points)
#endif
        for (int i=0; i<datasets[dataset]->NumOfPoints(); i++)
        {
            // Progress indication during fitting.
            #pragma omp critical
            cout << "\r" << i+1 << "/" << totalpoints << flush;

            // Index for this in the final data array
            int dataind=0;
            for (int dseti=0; dseti < dataset; dseti++)
                dataind += datasets[dseti]->NumOfPoints();
            dataind += i;
            
            
            
            double xbj      = datasets[dataset]->xbj(i);
            double y        = datasets[dataset]->y(i);              // inelasticity
            double Q2       = datasets[dataset]->Qsqr(i);
            double Q        = sqrt(Q2);
            double sigmar   = datasets[dataset]->ReducedCrossSection(i);
            double sigmar_err = datasets[dataset]->ReducedCrossSectionError(i);
            
            datavals[dataind] = sigmar;
            dataerrs[dataind] = sigmar_err;
            var_xbj[dataind] = xbj;
            var_y[dataind] = y;
            var_qsqr[dataind] = Q2;

            double theory=0, theory_charm=0;
            int calccount=0;
            if (!computeNLO && !useMasses) // Compute reduced cross section using leading order impact factors
            {
                theory = (fitsigma0)*SigmaComputer.SigmarLO(Q , xbj , y );
                ++calccount;
            }
            if (!computeNLO && useMasses)
            {
                theory=0;
                if (datasets[dataset]->OnlyCharm(i)==false)
                // if (true)
                {
                    theory = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj , y, false );
                }
                if (xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) < 0.01 and useCharm)
                {
                    theory_charm = (fitsigma0)*SigmaComputer.SigmarLOmass(Q , xbj*(1.0 + 4.0*1.35*1.35/(Q*Q)) , y, true ); // charm
                    if (datasets[dataset]->OnlyCharm(i) == true)
                        theory = theory_charm;
                    else
                        theory = theory + theory_charm;
                }
                  ++calccount;
            }
            
            if (computeNLO && !UseSub) // UNSUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useMasses){
                    if (useBoundLoop){
                        cout << "Z_2 bound dipole term not implemented yet with quark masses. EXIT." << endl;
                        exit(1);
                        // theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_UniformZ2Bound(Q , xbj , y );
                        ++calccount;}
                    if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                        // cout << "This should be in use???" << "charm mass is: " << qMass_charm << endl;
                        if (nlodis_config::MASS_SCHEME == nlodis_config::CHARM_ONLY){
                            theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_charm );
                        } else if (nlodis_config::MASS_SCHEME == nlodis_config::BEAUTY_ONLY){
                            theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_b_var );
                        } else if (nlodis_config::MASS_SCHEME == nlodis_config::LIGHT_PLUS_CHARM){
                            theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y )
                                     + (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_charm );
                        } else if (nlodis_config::MASS_SCHEME == nlodis_config::LIGHT_PLUS_CHARM_AND_BEAUTY){
                            theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y )
                                     + (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_charm )
                                     + (fitsigma0)*SigmaComputer.SigmarNLOunsub_massive(Q , xbj , y, qMass_b_var );
                        }
                        ++calccount;}
                }
                if (!useMasses){
                    if (useBoundLoop){
                        theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub_UniformZ2Bound(Q , xbj , y );
                        ++calccount;}
                    if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                        theory = (fitsigma0)*SigmaComputer.SigmarNLOunsub(Q , xbj , y );
                        //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
                        ++calccount;}
                    if (UseSigma3){
                        theory += (fitsigma0)*SigmaComputer.SigmarNLOunsub_sigma3(Q , xbj , y );
                        }
                }
            }
            
            if (computeNLO && UseSub) // SUB SCHEME Full NLO impact factors for reduced cross section
            {
                if (useMasses){
                    cout << "SUB SCHEME NOT IMPLEMENTED WITH QUARK MASSES. EXIT." << endl;
                    exit(1);
                }
                if (useBoundLoop){
                    theory = (fitsigma0)*SigmaComputer.SigmarNLOsub_UniformZ2Bound(Q , xbj , y );
                    ++calccount;}
                if (!useBoundLoop){ // the old way, no z2 lower bound in dipole loop term.
                    theory = (fitsigma0)*SigmaComputer.SigmarNLOsub(Q , xbj , y );
                    //theory = (fitsigma0)*SigmaComputer.SigmarNLOsubRisto(Q , xbj , y );
                    ++calccount;}
            }
            
            thdata[dataind] = theory;


            if (calccount>1)
            {
              cerr << "ERROR: Multiple computations. abort." << "count="<< calccount << endl;
              theory = 99999999;
              thdata[dataind] = 99999999;
              exit(1);
            }

            if (std::isnan(theory) or std::isinf(theory))
            {
                cerr << "Warning: theory result " << theory << " with parameters " << PrintVector(par) << endl;
                theory = 99999999;
            }

            //chisqr += datasets[dataset]->Weight()*SQR( (theory+theory_charm - sigmar) / sigmar_err );
            points = points + datasets[dataset]->Weight();
        }
    }
    cout << endl;
    // Minimize sigma02
    std::vector<double> sigma02fit = FindOptimalSigma02(datavals,dataerrs, thdata);
    double chisqr_over_n = sigma02fit[1];
    double sigma02 = sigma02fit[0];
    // comparing to old sigma02
    double old_chisq_over_n = minimiser_helper_chisqr_vec(old_sigma02, datavals, dataerrs, thdata);
    // Output for plotting
    cout << " ## var_xbj, var_qsqr, var_y, sigmar, sigmar_err, sigma02xthdata " << endl;
    if(nlodis_config::PRINTDATA){
        for(int i=0; i<var_xbj.size(); i++){
        #pragma omp critical
        cout    << setw(10) << var_xbj[i]      << " "
                << setw(10) << var_qsqr[i]     << " "
                << setw(10) << var_y[i]        << " "
                << setw(10) << datavals[i]     << " "
                << setw(10) << dataerrs[i]     << " "
                << setw(10) << old_sigma02*thdata[i] /* << " "
                << setw(10) << sigma02*theory_charm <<*/ 
                << endl;
            }
    }
    cout << " ## sigma r data ends here " << endl;    
    // cout    << endl 
    //         << "# Calculated chi^2/N = " << chisqr_over_n
    //         << " (N=" << points
    //         << "), parameters (" << PrintVector(par)
    //         << ", sigma02=" << sigma02
    //         << ")" << endl<<endl;
    cout    << endl
            << "# Calculated chi^2/N = " << old_chisq_over_n
            << " (N=" << points
            << "), parameters (" << PrintVector(par)
            << ", sigma02=" << old_sigma02
            << ")" << endl<<endl;
    

    // delete DipolePointer;
    return chisqr_over_n * points;  // Return chisqr for Minuit error estimation.
}

void NLODISFitter::AddDataset(Data& d)
{
    datasets.push_back(&d);
}

NLODISFitter::NLODISFitter(MnUserParameters parameters_)
{
    parameters = parameters_;
    
    cubaMethod = "suave";  // Default choise for Cuba
}

string PrintVector(vector<double> v)
{
    stringstream ss;
    for (int i=0; i <v.size(); i++)
    {
        ss << v[i];
        if (i < v.size()-1)
            ss <<", ";
    }
    return ss.str();
}

///===========================================================================================
// CUBA WRAP
namespace cuba_config{
    int verbose=0;
    //maxeval=nlodis_config::CUBA_MAXEVAL;
    //double epsrel=nlodis_config::CUBA_EPSREL;
    double epsabs=0;
}

void Cuba(string method, int ndim, integrand_t integrand,
    void *userdata, double *integral, double *error, double *prob) {
    // common arguments
    int ncomp=1, nvec=1, seed=0, mineval=0, last=4;
    int nregions, neval, fail;
    void *spin=NULL;
    char *statefile=NULL;
    if(method=="vegas"){
    // Vegas-specific arguments
    int nstart=1000, nincrease=500, nbatch=1000, gridno=0;
    Vegas(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose,seed,mineval,
        nlodis_config::CUBA_MAXEVAL,nstart,nincrease,nbatch,gridno,statefile,
        spin,&neval,&fail,integral,error,prob);
    }
    else if(method=="suave"){
    // Suave-specific arguments
    int nnew=1e3, nmin=2; // nnew=10e3
    double flatness=25; //25;
    Suave(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose | last,seed,mineval,
        nlodis_config::CUBA_MAXEVAL,nnew,nmin,flatness,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="divonne"){
    if(ndim==1) ndim=2;
    // Divonne-specific arguments
    int key1=1*47, key2=1, key3=1, maxpass=5, ngiven=0, nextra=0;
    double border=1e-8, maxchisq=10, mindeviation=0.25;
    Divonne(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose,seed,mineval,
        nlodis_config::CUBA_MAXEVAL,key1,key2,key3,maxpass,border,maxchisq,
        mindeviation,ngiven,ndim,NULL,nextra,NULL,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
    else if(method=="cuhre"){
    if(ndim==1) ndim=2;
    // Cuhre-specific arguments
    int key=0;
    Cuhre(ndim,ncomp,integrand,userdata,nvec,nlodis_config::CUBA_EPSREL,
        cuba_config::epsabs,cuba_config::verbose | last,mineval,
        nlodis_config::CUBA_MAXEVAL,key,statefile,spin,
        &nregions,&neval,&fail,integral,error,prob);
    }
}




///===========================================================================================
///===========================================================================================
// HELPERS
double ComputeSigmaR::Sr(double r, double x) {
    double Srx;//, Nrx;
    // cout << "Sr r x: " << r << " " << x << endl;
    if(r<nlodis_config::MINR){
        Srx = 1.;
    }else if(r>nlodis_config::MAXR-1e-7){
        Srx = 0.;
    }else{
        // Srx = ClassScopeDipolePointer->S_y(r, log(1/x))
        // Note: icY0 controls how much evolution we have before we define that we have "initial condition"
        if (x > icX0_bk){
            Srx = ClassScopeDipolePointer->S(r, (icX0_bk*std::exp(-icY0)) ) ; //1-Nrx;
        } else {
            Srx = ClassScopeDipolePointer->S(r, (x*std::exp(-icY0)) ) ; //1-Nrx;
        }
    }
    return Srx;
}

double ComputeSigmaR::SrY(double r, double Y) {
    double Sry;
    if(r<nlodis_config::MINR){
        Sry = 1.;
    }else if(r>nlodis_config::MAXR-1e-7){
        Sry = 0.;
    }else{
        // Note: icY0 controls how much evolution we have before we define that we have "initial condition"
        if (std::exp(-Y) > icX0_bk){
            Sry = ClassScopeDipolePointer->S(r, (icX0_bk*std::exp(-icY0)) ) ; //1-Nrx;
        } else {
            Sry = ClassScopeDipolePointer->S(r, (std::exp(-Y)*std::exp(-icY0)) ) ; //1-Nrx;
        }
    }
    return Sry;
}

double ComputeSigmaR::SrTripole(double x01, double x_x01, double x02, double x_x02, double x21, double x_x21) {
    // double x_x01 = cov_to_eta_x(x01, x);
    // double x_x02 = cov_to_eta_x(x02, x);
    // double x_x21 = cov_to_eta_x(x21, x);
    return Nc/(2*CF)*(Sr(x02, x_x02)*Sr(x21, x_x21) - 1/Sq(Nc)*Sr(x01, x_x01));
}

double ComputeSigmaR::P(double z) {
    return 1.0-z+Sq(z)*0.5; // EI TÄSMÄÄ ARTIKKELIN P-funktion KANSSA. P=1/2*Partikkeli
}

double ComputeSigmaR::heaviside_theta(double x) {
    double result = 0;
    if (x >= 0) {result = 1.0;}
    return result;
}



///===========================================================================================
// LO Rapidities with / w/o shift
double ComputeSigmaR::Xrpdty_LO_projectileY(double x, double Qsq, double rsq){
    double X = (this->*Xrpdty_LO_projectileY_z2min_PTR)( x, Qsq );
    return X;
}

double ComputeSigmaR::Xrpdty_LO_targetETA(double x, double Qsq, double rsq){
    // Target rapidity eta = Y - rho is shifted so here we compute the bjorken-x corresponding to the shifted rapidity.
    // x_eta = x_Y * exp(rho)
    // THIS SHOULD ONLY be used with sub scheme LO term and with the DIPOLE term. Not with unsub IC.

    double X_Y = (this->*Xrpdty_LO_projectileY_z2min_PTR)( x, Qsq );
    // double rho = rho_rapidity_shift_QQ0(Qsq);
    double rho = rho_rapidity_shift(Qsq, rsq);
    return X_Y*std::exp(rho); // return X_eta for the shifted eta = Y - rho
}


// NLO Rapidities and shift(s)
double ComputeSigmaR::Xrpdty_NLO_projectileY(double Qsq, double z2, double z2min, double icX0, double x01sq = 0, double x02sq = 0, double x21sq = 0 ){
    double X = z2min * icX0/z2;
    return X;
}
double ComputeSigmaR::Xrpdty_NLO_targetETA(double Qsq, double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq ){
    // Target rapidity eta = Y - rho is shifted so here we compute the bjorken-x corresponding to the shifted rapidity.
    // x_eta = x_Y * exp(rho)
    double X_Y = Xrpdty_NLO_projectileY(Qsq, z2, z2min, icX0);
    // cout << "This QQG term level uniform shift is wrong! Don't use this!" << endl;
    // TODO test that this is used only with QQ0 shift
    // double rho = rho_rapidity_shift(Qsq, x01sq, x02sq, x21sq);
    // return X_Y*std::exp(rho);
    return X_Y;
}

double ComputeSigmaR::cov_to_eta_x(double r, double xbj, double Qsq){
    // shift xbj_Y to xbj_eta
    double rho = rho_rapidity_shift(Qsq, Sq(r));
    return xbj * std::exp(rho);
}

// Target rapidity Eta shift calculator(s)
double ComputeSigmaR::rho_no_shift_projectileY(double Qsq, double x01sq, double x02sq, double x21sq){
    return 0;
}

double ComputeSigmaR::rho_rapidity_shift_QQ0(double Qsq, double x01sq, double x02sq, double x21sq){
    double rho;
    rho = std::log( Qsq/icQ0sqr);

    if ( rho < 0 ){
        return 0;
    } else {
        return rho;
    }
}

double ComputeSigmaR::rho_rapidity_shift_RQ0(double Qsq, double Rsq, double x02sq, double x21sq){
    // This is the dipole size dependent target eta shift; Rsq is the relevant dipole: parent or daughter.
    double rho;
    rho = std::log( 1/(Rsq * icQ0sqr));

    if ( rho < 0 ){
        return 0;
    } else {
        return rho;
    }
}

double ComputeSigmaR::rho_rapidity_shift_XR(double Qsq, double x01sq, double x02sq, double x21sq){
    double rho;
    rho = std::log( 1/icQ0sqr * x02sq/x01sq);

    if ( rho < 0 ){
        return 0;
    } else {
        return rho;
    }
}

double ComputeSigmaR::rho_rapidity_shift_MAX_XYR(double Qsq, double x01sq, double x02sq, double x21sq){
    double rho;
    rho = std::log( 1/icQ0sqr * std::max(x02sq, x21sq)/x01sq);

    if ( rho < 0 ){
        return 0;
    } else {
        return rho;
    }
}

double ComputeSigmaR::x_eta_delta_ij_r(double x_ij_sq, double rsq){
    // Edmond et al 1902.06637 Eq. 5.7
    // Since this is used as a negative shift eta' = (eta - delta_ij) one gets the corresponding shifted x_eta' as x_eta/x_delta_ij.
    double delta_ij_r = std::max(0., std::log(rsq / x_ij_sq));
    double x_eta_delta = std::exp(-delta_ij_r);
    return x_eta_delta;
}


///===========================================================================================
// RUNNING COUPLINGS
double ComputeSigmaR::alpha_bar_running_pd_sharpcutoff( double rsq ) { // alphabar = Nc/M_PI * alphas
    double logvar = std::log(	(4*(alpha_scaling_C2_)/(rsq*Sq(lambdaqcd))	)	);
    double AlphaSres = 12*M_PI/(	(33.0-2*Nc)*logvar );
    if ((logvar<0) || (AlphaSres>0.7)) {
        AlphaSres = 0.7;
        }
    return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_running_pd( double rsq ) { // alphabar = Nc/M_PI * alphas
    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(rsq*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
    return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_fixed( double rsq ) { // alphabar = Nc/M_PI * alphas
    return 0.190986; // alphas_fixed=0.2
}


// QG term specific running couplings
struct Alphasdata{
    double x01sq,x02sq,x21sq;
};

double ComputeSigmaR::alpha_bar_QG_fixed( void *userdata ) { // alphabar = Nc/M_PI * alphas
    return 0.190986; // alphas_fixed=0.2
}

double ComputeSigmaR::alpha_bar_QG_running_pd( void *userdata ) { // alphabar = Nc/M_PI * alphas
    Alphasdata *dataptr = (Alphasdata*)userdata;
    double rsq=dataptr->x01sq;
    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(rsq*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
  return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_QG_running_guillaume( void *userdata ) { // alphabar = Nc/M_PI * alphas
    Alphasdata *dataptr = (Alphasdata*)userdata;
    double x01sq=dataptr->x01sq;
    double x02sq=dataptr->x02sq;
    double x21sq=dataptr->x21sq;
    double r_eff_sqr = x01sq * std::pow( x21sq / (x02sq), (x02sq-x21sq)/(x01sq) ); // r = x01 , X = x02 , Y = x21 , Q_{123} = 4C^2 / r_eff_sqr.

    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(r_eff_sqr*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
  return Nc/M_PI*AlphaSres;
}

double ComputeSigmaR::alpha_bar_QG_running_smallest( void *userdata ) { // alphabar = Nc/M_PI * alphas
    Alphasdata *dataptr = (Alphasdata*)userdata;
    double x01sq=dataptr->x01sq;
    double x02sq=dataptr->x02sq;
    double x21sq=dataptr->x21sq;
    double r_min_0 = std::min(x01sq, x02sq); // r = x01 , X = x02 , Y = x21 , Q_{123} = 4C^2 / r_eff_sqr.
    double r_min_sqr = std::min(r_min_0, x21sq);

    double scalefactor = 4.0*alpha_scaling_C2_;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;

    double AlphaSres = 4.0*M_PI / ( b0 * std::log(
    std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(scalefactor/(r_min_sqr*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
    ) );
  return Nc/M_PI*AlphaSres;
}


///===========================================================================================
// HIGHER ORDER BK SIGMA3 CORRECTION TERM KERNELS

// resum bk
double ComputeSigmaR::K_resum (double rsq, double x20sq, double x21sq) {
    // returns K_DLA * K_STL
    // double log resummations r = x01 , X = x02 , Y = x21
    double xsq =  std::log(x20sq/rsq) * std::log(x21sq/rsq) ;
    double K_DLA;
    double resummation_alphabar = Alphabar(rsq);
    if (xsq >=0)
        {            
            // argument to the Bessel function is 2sqrt(bar as * xsq)
            double as_x = std::sqrt( resummation_alphabar * xsq );
            K_DLA = gsl_sf_bessel_J1(2.0*as_x) / as_x;
        }
    else // L_xzr L_yzr < 0
        {
            xsq = std::abs(xsq);
            double as_x = std::sqrt( resummation_alphabar * xsq);
            gsl_sf_result res;
            int status = gsl_sf_bessel_I1_e(2.0*as_x, &res);
            
            if (status != GSL_SUCCESS)
            {
                if (std::isnan(x20sq) or std::isnan(x21sq))
                    return 0;	// 0/0, probably z=x or z=y
                cerr << "GSL error " << status <<", result " << res.val << ", as_x=" << as_x << ", xsq=" << xsq << LINEINFO << endl;
                return 0;
            }
            K_DLA = res.val / as_x; 
        }
    if (std::isnan(K_DLA))
        {
            K_DLA =  1;
        }
    // single log resummations
    double K_STL;
    double minxysq = std::min(x20sq,x21sq);
    if (minxysq < 1e-40) minxysq = 1e-40;

    const double A1 = 11.0/12.0;
    K_STL = std::exp( - resummation_alphabar * A1 * std::abs( std::log( config::KSUB * rsq/minxysq ) ) );

    double K_RES = K_DLA*K_STL; // 1601.06598
    return K_RES;
}

double ComputeSigmaR::K_lobk (double rsq, double x20sq, double x21sq) {
    double K_loBK = rsq/(x20sq*x21sq); // 1601.06598
    return K_loBK;
}



///===========================================================================================
// NLO sub scheme qg subtraction terms for different evolution equations.
double ComputeSigmaR::ILNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors numerically
    double subterm;
    subterm = ILNLOqg(Q,x,z1, 0,x01sq,x02sq,x21sq);
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors explicitly by hand
    double subterm;
    double impact_factor_lo = ILLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    double x01 = sqrt(x01sq);
    double x02 = sqrt(x02sq);
    double x21 = sqrt(x21sq);
    double nlo_if_bk_evol_dipoles = Sr(x01, x) - SrTripole(x01, x, x02, x, x21, x);
    subterm = impact_factor_lo * nlo_if_kernel * nlo_if_bk_evol_dipoles;
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the RESUM kernel
    double subterm;
    double K_res = K_resum(x01sq, x02sq, x21sq);
    double lobk_subterm = ILNLOqg_subterm_lobk_explicit(Q,x,z1,z2,x01sq,x02sq,x21sq);
    subterm = K_res * lobk_subterm;
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the KCBK kernel from 1401.0313
    double delta012 = std::max(0.0, std::log( std::min(x02sq, x21sq) / (x01sq) ) ); // (166)
    // double shifted_rapidity = helper->rapidity - delta012;
    double scaled_xbj = x * std::exp(delta012); // rapidity shift Y' = Y - Delta012 leads to a scaling of xbj
    if (scaled_xbj > icX0) // Heaviside step with Y_shift > 0 leads to xbj_scaled < 0 and so the cutoff inequality flips.
        return 0;   // Step function in (165)
    
    // Dipoles at shifted rapidity / scaled xbj
    double s02_shifted = Sr(sqrt(x02sq), scaled_xbj);
    double s12_shifted = Sr(sqrt(x21sq), scaled_xbj);
    double s01 = Sr(sqrt(x01sq), x);
    
    double subterm;
    double impact_factor_lo = ILLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    subterm = impact_factor_lo * nlo_if_kernel * ( -s02_shifted * s12_shifted + s01);
    return subterm;
}

double ComputeSigmaR::ILNLOqg_subterm_trbk_edmond(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // subtraction term for the Edmond et al. kinematically constrained target rapidity BK. arxiv 1902.06637 Eq. 9.3.
    // N.B. No step function in the final form Eq. 9.3.; Only LOBK kernel and the dipole amplitudes at shifted rapidities.
    // Thus we only need to add the LO impact factor here then.

    // Shifted rapidites; delta defined in Eq. 5.7
    double z2min = z2lower_bound(x,Q*Q);
    double x_eta = Xrpdty_NLO_targetETA(Q*Q, z2, z2min, icX0, x01sq, x02sq, x21sq);
    double x_eta_delta02r = x_eta / x_eta_delta_ij_r(x02sq, x01sq);
    double x_eta_delta21r = x_eta / x_eta_delta_ij_r(x21sq, x01sq);
    
    // Dipoles at shifted rapidity / scaled xbj
    double s02_shifted = Sr(sqrt(x02sq), x_eta_delta02r);
    double s12_shifted = Sr(sqrt(x21sq), x_eta_delta21r);
    double s01 = Sr(sqrt(x01sq), x_eta);

    double subterm;
    double impact_factor_lo = ILLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    subterm = impact_factor_lo * nlo_if_kernel * ( -s02_shifted * s12_shifted + s01);
    return subterm;
}



double ComputeSigmaR::ITNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors numerically
    double subterm;
    subterm = ITNLOqg(Q,x,z1,0,x01sq,x02sq,x21sq);
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit of the full NLO impact factors explicitly by hand
    double subterm;
    double impact_factor_lo = ITLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq);
    double x01 = sqrt(x01sq);
    double x02 = sqrt(x02sq);
    double x21 = sqrt(x21sq);
    double nlo_if_bk_evol_dipoles = Sr(x01, x) - SrTripole(x01, x, x02, x, x21, x);
    subterm = impact_factor_lo * nlo_if_kernel * nlo_if_bk_evol_dipoles;
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the RESUM kernel
    double subterm;
    double K_res = K_resum(x01sq, x02sq, x21sq);
    double lobk_subterm = ITNLOqg_subterm_lobk_explicit(Q,x,z1,z2,x01sq,x02sq,x21sq);
    subterm = K_res * lobk_subterm;
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // z2 -> 0 limit LOBK kernel complemented with the KCBK kernel from 1401.0313
    double delta012 = std::max(0.0, std::log( std::min(x02sq, x21sq) / (x01sq) ) ); // (166)
    // double shifted_rapidity = helper->rapidity - delta012;
    double scaled_xbj = x * std::exp(delta012); // rapidity shift Y' = Y - Delta012 leads to a scaling of xbj
    if (scaled_xbj > icX0) // Heaviside step with Y_shift > 0 leads to xbj_scaled < 0 and so the cutoff inequality flips.
        return 0;   // Step function in (165)
    
    // Dipoles at shifted rapidity / scaled xbj
    double s02_shifted = Sr(sqrt(x02sq), scaled_xbj);
    double s12_shifted = Sr(sqrt(x21sq), scaled_xbj);
    double s01 = Sr(sqrt(x01sq), x);
    
    double subterm;
    double impact_factor_lo = ITLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    subterm = impact_factor_lo * nlo_if_kernel * ( -s02_shifted * s12_shifted + s01);
    return subterm;
}

double ComputeSigmaR::ITNLOqg_subterm_trbk_edmond(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
    // subtraction term for the Edmond et al. kinematically constrained target rapidity BK. arxiv 1902.06637 Eq. 9.3.
    // N.B. No step function in the final form Eq. 9.3.; Only LOBK kernel and the dipole amplitudes at shifted rapidities.
    // Thus we only need to add the LO impact factor here then.

    // Shifted rapidites; delta defined in Eq. 5.7
    double z2min = z2lower_bound(x,Q*Q);
    double x_eta = Xrpdty_NLO_targetETA(Q*Q, z2, z2min, icX0, x01sq, x02sq, x21sq);
    double x_eta_delta02r = x_eta / x_eta_delta_ij_r(x02sq, x01sq);
    double x_eta_delta21r = x_eta / x_eta_delta_ij_r(x21sq, x01sq);
    
    // Dipoles at shifted rapidity / scaled xbj
    double s02_shifted = Sr(sqrt(x02sq), x_eta_delta02r);
    double s12_shifted = Sr(sqrt(x21sq), x_eta_delta21r);
    double s01 = Sr(sqrt(x01sq), x_eta);

    double subterm;
    double impact_factor_lo = ITLO(Q,z1,x01sq);
    double P_at_zero = 2;
    double nlo_if_kernel = P_at_zero * 0.5 * (x01sq + x21sq - x02sq) / (x02sq * x21sq); // this is the symmetric half of the LOBK kernel
    subterm = impact_factor_lo * nlo_if_kernel * ( -s02_shifted * s12_shifted + s01);
    return subterm;
}



///===========================================================================================
// Computation dipole & data passing helpers
ComputeSigmaR::ComputeSigmaR(AmplitudeLib *ObjectPointer){
    ClassScopeDipolePointer = ObjectPointer;
}

struct Userdata{
    ComputeSigmaR* ComputerPtr;
    double Q, xbj;
    double qMass;
    double icX0;
};


double sumef_from_mass(double mf){
    double sumef; // \sum_f e_f^2
    // return 1.0; // mass function plots
    if (mf < 1e-8) {
        sumef = 6.0/9.0; // light quarks uds only.
    } else if (mf > 0.7*0.00216 && mf < 1.3*0.00216) {
        // u
        sumef = 4.0/9.0;
    } else if (mf > 0.7*0.00467 && mf < 1.3*0.00467) {
        // d
        sumef = 1.0/9.0;
    } else if (mf > 0.7*0.093 && mf < 1.3*0.093) {
        // s
        sumef = 1.0/9.0;
    } else if (mf > 0.7*1.27 && mf < 2*1.27) {
        // c
        sumef = 4.0/9.0;
    } else if (mf > 0.7*4.180 && mf < 3*4.180) {
        // b
        sumef = 1.0/9.0;
    } else {
        // failure to identify quark from mass: kill it all.
        cout << "Failure to identify quark sumef from mass. Exit." << endl;
        exit(1);
    }
    return sumef;
}

///===========================================================================================
///===========================================================================================
///===========================================================================================
// INTEGRATORS for cross sections

/*
// --- L L L --- LO -------- L L L --- LO -------- L L L --- LO -----------
*/
double ComputeSigmaR::ILLO(double Q, double z1, double x01sq) {
    double bessel_inner_fun = Q*sqrt(z1*(1.0-z1)*x01sq);
    double res = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res = 0;
    }else{
        res = 4.0*Sq(Q)*Sq(z1)*Sq(1.0-z1)*Sq(gsl_sf_bessel_K0( bessel_inner_fun ));
    }   
    return res;
}

// MASSLESS
int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=x01*x01;
    double Xrpdty_lo = Optr->Xrpdty_LO(xbj, Sq(Q), x01sq);

    double res=(1.0-(Optr->Sr(x01,Xrpdty_lo)))*(Optr->ILLO(Q,z1,x01sq))*x01;
        *f=res;
    return 0;
}

double ComputeSigmaR::LLOp(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Cuba(cubamethod,ndim,integrand_ILLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

// MASS
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double qmass=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double Xrpdty_lo = Optr->Xrpdty_LO(xbj, Sq(Q), x01*x01);

    double af = sqrt( Sq(Q)*z1*(1-z1) + Sq(qmass) );
    double impactfac = 4.0*Sq(Q*(z1)*(1.0-z1)*gsl_sf_bessel_K0(af*x01));
    double res=(1.0-(Optr->Sr(x01,Xrpdty_lo)))*(impactfac)*x01;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::LLOpMass(double Q, double x, bool charm) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    if (charm==false)
        userdata.qMass=qMass_light;
    else
        userdata.qMass=qMass_charm;
    userdata.ComputerPtr=this;
    double ef = sumef;
    if (charm) ef = 4.0/9.0;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*ef;
    Cuba(cubamethod,ndim,integrand_ILLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LLOp_massive(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;

    double ef = sumef_from_mass(mf);

    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*ef;
    Cuba(cubamethod,ndim,integrand_ILLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

///===========================================================================================
/*
// --- T T T --- LO -------- T T T --- LO -------- T T T --- LO -----------
*/

double ComputeSigmaR::ITLO(double Q, double z1, double x01sq) {
    double bessel_inner_fun = Q*sqrt(z1*(1.0-z1)*x01sq);
    double res = 0;
    if (bessel_inner_fun < 1e-7){
        // cout << bessel_inner_fun << " " << Q << " " << z1 << " " << x01sq << endl;
        res = 0;
    }else{
        res = Sq(Q)*z1*(1.0-z1)*(1.0-2.0*z1+2.0*Sq(z1))*Sq(gsl_sf_bessel_K1(bessel_inner_fun));
    }
    return res;
}

// MASSLESS
int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=x01*x01;
    double Xrpdty_lo = Optr->Xrpdty_LO(xbj, Sq(Q), x01sq);

    double res=(1.0-(Optr->Sr(x01,Xrpdty_lo)))*(Optr->ITLO(Q,z1,x01sq))*x01;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::TLOp(double Q, double x) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Cuba(cubamethod,ndim,integrand_ITLOp,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

// MASS
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double qmass=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double Xrpdty_lo = Optr->Xrpdty_LO(xbj, Sq(Q), x01*x01);

    double af = sqrt( Sq(Q)*z1*(1.0-z1) + Sq(qmass) );
    double impactfac = (1.0-2.0*z1+2.0*Sq(z1))*Sq(af*gsl_sf_bessel_K1(af*x01)) + Sq( qmass*gsl_sf_bessel_K0( af*x01 ) );
    double res=(1.0-(Optr->Sr(x01,Xrpdty_lo)))*(impactfac)*x01;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::TLOpMass(double Q, double x, bool charm) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    if (charm==false)
        userdata.qMass=qMass_light;
    else
        userdata.qMass = qMass_charm;
    userdata.ComputerPtr=this;

    double ef = sumef;
    if (charm) ef = 4.0/9.0;

    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*ef;
    Cuba(cubamethod,ndim,integrand_ITLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TLOp_massive(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=2;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass = mf;
    userdata.ComputerPtr=this;

    double ef = sumef_from_mass(mf);

    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*ef;
    Cuba(cubamethod,ndim,integrand_ITLOpMass,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

///===========================================================================================
///===========================================================================================
///===========================================================================================
/*
// --- L L L --- NLO -------- L L L --- NLO -------- L L L --- NLO -----------
*/

///*
// NLO
double ComputeSigmaR::ILdip(double Q, double z1, double x01sq) { //ILbeufDIP
    double fac1 = Sq(z1)*Sq(1.0 - z1);
    double facNLO = 4.0*Sq(Q*gsl_sf_bessel_K0(Q*sqrt(z1*(1-z1)*x01sq)));
    return fac1*facNLO;
}

int integrand_ILdip(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) { //integrand_ILbeufDIP
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->Alphabar(x01sq); //2;
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double regconst = 5.0/2.0 - Sq(M_PI)/6.0;
    double res;

    res = SKernel*(Optr->ILdip(Q,z1,x01sq))*x01*( alphfac*(0.5*Sq(log(z1/(1-z1))) + regconst ) );
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILdip_z2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    double res;
    if (z2min > 0.5){ // Check that z2min is not too large since z1 limits are [z2min , 1-z2min]. IF it is too large, return res=0.
        res = 0;
        *f=res;
        return 0;
    }
    double z1=((1.0-z2min)-z2min)*x[0]+z2min;
    double z2=((z1)-z2min)*x[1]+z2min; // z2 integration variable for k1-term AND for the mirror symmetry combined integral!
    double x01=nlodis_config::MAXR*x[2];
    double x01sq=Sq(x01);
    double jac=((1.0-z2min)-z2min)*((z1)-z2min); // jacobiaaniin skaalaus molemmista z-muuttujista!
    double Xrpdt= z2min * X0/z2; // consistent with the qg-terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdt); //ok?

    double fivefourths = 5.0/4.0;
    double loopcontribution = 1.0/z2 * ( log(1.0 + z2 / (1.0-z1) ) + log(1.0 - z2 / z1) ) + 1.0/z1 * fivefourths ; // a.k.a k1- or z1-term
    res = 2 * jac * SKernel*(Optr->ILdip(Q,z1,x01sq))*x01*( alphfac*( loopcontribution ) ); // factor of 2 is from the combination of the mirror symmteric q/qbar terms.

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double x01=sqrt(x01sq);
    double x02=sqrt(x02sq);
    double x21=sqrt(x21sq);

    double X3sq = z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;

    double Qsq = Sq(Q);
    double x_x01 = cov_to_eta_x(x01, x, Qsq);
    double x_x02 = cov_to_eta_x(x02, x, Qsq);
    double x_x21 = cov_to_eta_x(x21, x, Qsq);

    double bessel_innerfun = Q*sqrt(X3sq);
    double facNLO = 0;
    // if (bessel_innerfun < 1e-7){
    if (bessel_innerfun > 200){
        // cout << "bessel_innerfun = " << bessel_innerfun << " Q " << Q << " X3sq " << X3sq << endl;
        return 0;
    } else if (bessel_innerfun < 1e-30){
        return 1;
    }else{
        facNLO = 4.0*Sq(Q*gsl_sf_bessel_K0( bessel_innerfun ));
    }
    return facNLO*(1-SrTripole(x01, x_x01, x02, x_x02, x21, x_x21));
}

double ComputeSigmaR::ILNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) { // old ILNLObeufQG
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);

    double fac1 = Sq(z1)*Sq(1.0 - z1);
    double xi 	= z2/(1.0-z1);
    double fun1 = 1+Sq(1-xi);
    double fun2 = x20x21/(x02sq*x21sq);
    double fac2 = 1/x02sq - fun2;
    double fac3 = Sq(xi)*fun2;

    double facNLO1 = Bessel0Tripole(Q, x, z1, z2, x01sq, x02sq, x21sq);
    double facNLO2 = Bessel0Tripole(Q, x, z1, z2, x01sq, 0    , x01sq);

    double res = fac1*((fun1*fac2)*(facNLO1 - facNLO2) + fac3*facNLO1 );
    if(gsl_finite(res)==1){
        return res;
    }else{
        res=0;
    }
    return res;
}

double ComputeSigmaR::ILNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double impactfac_lo = ILLO(Q,z1,x01sq);
    double BKkernel = (K_kernel(x01sq,x02sq,x21sq) - 1.0)*2*K_lobk(x01sq,x02sq,x21sq);
    double S012 = SrTripole(sqrt(x01sq),x,sqrt(x02sq),x,sqrt(x21sq),x);
    double Sdipole_bkevol = Sr(sqrt(x01sq),x) - S012;

    double res = BKkernel*Sdipole_bkevol*impactfac_lo;
    return res;
}

double ComputeSigmaR::ILNLOqgRisto(double Q, double x, double z0, double z2, double x01sq, double x02sq, double x21sq){ // old ILNLOristo2QG
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double z3 = 1 - z0 - z2;
    double X3sq = z0*z3*x01sq + z0*z2*x02sq + z2*z3*x21sq;

    double fun2 = x20x21/(x02sq*x21sq);

    // FIRST FORM
    double facterm1 = Sq(z3) * ( 2*z0*(z0+z2) + Sq(z2) )/x02sq;
    double facterm2 = Sq(z0) * ( 2*z3*(z2+z3) + Sq(z2) )/x21sq;
    double facterm3 = 2*( (z0+z2)*z0*Sq(z3) + (z2+z3)*z3*Sq(z0) ) * fun2;

    /* // SECOND FORM
    double facterm1 = Sq(z3) * ( 2*z0*(z0+z2) + Sq(z2) ) * ( 1/(x02sq) - fun2 );
    double facterm2 = Sq(z0) * ( 2*z3*(z2+z3) + Sq(z2) ) * ( 1/(x21sq) - fun2 );
    double facterm1expn = Sq(z3) * ( 2*z0*(z0+z2) + Sq(z2) )*(1/x02sq);
    double facterm2expo = Sq(z0) * ( 2*z3*(z2+z3) + Sq(z2) )*(1/x21sq);
    double facterm3 = Sq(z2)*( Sq(z0) + Sq(z3) ) * fun2;
    */

    // double facBesselTrip = Bessel0Tripole(Q, x, z3, z2, x01sq, x02sq, x21sq, X3sq);
    double facBesselTrip = Bessel0Tripole(Q, x, z3, z2, x01sq, x02sq, x21sq);

    double Qbarn = sqrt(z3*(1-z3))*Q;
    double Qbaro = sqrt(z0*(1-z0))*Q;
    double SKernel = 1.0 - Sr(sqrt(x01sq),x);

    ///*
    double facExpn, facExpo;
    if(x02sq/x01sq>1e-8 && x02sq/x01sq<5e2){
        facExpn = std::exp(-x02sq/(std::exp(M_EULER)*x01sq));
    }else if (x02sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }
    if(x21sq/x01sq>1e-8 && x21sq/x01sq<5e2){
        facExpo = std::exp(-x21sq/(std::exp(M_EULER)*x01sq));
    }else if (x21sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }

    double facBesExpn, facBesExpo;
    if(Qbarn*sqrt(x01sq)>1e-8 && Qbarn*sqrt(x01sq)<5e2){
        facBesExpn = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbarn*sqrt(x01sq)))*facExpn*SKernel;
    }else{ facBesExpn = 0; }
    if(Qbaro*sqrt(x01sq)>1e-8 && Qbaro*sqrt(x01sq)<5e2){
        facBesExpo = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbaro*sqrt(x01sq)))*facExpo*SKernel;
    }else{ facBesExpo = 0; }
    //*/

    //double facBesExpn = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbarn*sqrt(x01sq)))*std::exp(-x02sq/(std::exp(M_EULER)*x01sq))*SKernel;
    //double facBesExpo = 4.0*Sq(Q*gsl_sf_bessel_K0(Qbaro*sqrt(x01sq)))*std::exp(-x21sq/(std::exp(M_EULER)*x01sq))*SKernel;

    // FIRST FORM
    double term1 = facterm1*( facBesselTrip - facBesExpn );
    double term2 = facterm2*( facBesselTrip - facBesExpo );
    double term3 = facterm3*facBesselTrip;
    double res = term1 + term2 - term3;

    /* //SECOND FORM
    double term1 = facterm1* facBesselTrip - facterm1expn*facBesExpn;
    double term2 = facterm2* facBesselTrip - facterm2expo*facBesExpo;
    double term3 = facterm3*facBesselTrip;
    double res = term1 + term2 + term3;
    */

    return res;
}

int integrand_ILqgunsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) { //integrand_ILbeufQGiancu
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min; // linear z_2 int
    // double ymax = std::log(1-z1);
    // double ymin = std::log(z2min);
    // double y=(ymax-ymin)*x[1]+ymin;
    // double z2=std::exp(y);
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    // double rhomax=std::log(nlodis_config::MAXR);
    // double rhomin=std::log(nlodis_config::MINR);
    // double rho01=(rhomax-rhomin)*x[2] + rhomin;
    // double rho02=(rhomax-rhomin)*x[3] + rhomin;
    // double x01=std::exp(rho01);
    // double x02=std::exp(rho02);
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min)*Sq(nlodis_config::MAXR); // linear z_2 int, linear x01 x02
    // double jac=(1.0-z2min)*(ymax-ymin)*z2; // log z_2 int
    // double jac=(1.0-z2min)*(1.0-z1-z2min)*Sq(rhomax-rhomin)*x01*x02; // linear z_2 int, log x01, log x02

    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res =   jac*alphfac*( Optr->ILNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res = jac*alphfac*( Optr->ILNLOsigma3(Q,Xrpdt,z1,0,x01sq,x02sq,x21sq) )/z2*x01*x02; // K_NLO(z_2 = 0)

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min >= 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=x[0];
    double z2=(1.0-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double ILNLOqg_z2 = 0;
    if (1-z1-z2 > 0){ // heaviside theta front factor on the full NLO term
        ILNLOqg_z2 = (Optr->ILNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq));
    }
    double ILNLOqg_z2_to_0 = Optr->ILNLOqg_subterm(Q,Xrpdt,z1, 0,x01sq,x02sq,x21sq);

    double res = jac*alphfac*(
                    ILNLOqg_z2
                    -
                    ILNLOqg_z2_to_0
                )/(z2)*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=(1.0-x0lim)*x[0];
    double z2=(1.0-z0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim)*(1.0-z0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq); //0.2;
    double alphfac=alphabar*CF/Nc;

    double res;

    res =   jac*alphfac*( Optr->ILNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=x[0];
    double z2=(1.0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq); //0.2;
    double alphfac=alphabar*CF/Nc;

    double res;

    res =   jac*alphfac*(
                    (Optr->ILNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq))*(Optr->heaviside_theta(1-z0-z2))
                    -Optr->ILNLOqgRisto(Q,Xrpdt,z0, 0,x01sq,x02sq,x21sq)
                    )/z2*x01*x02;

        if(gsl_finite(res)==1){
            *f=res;
        }else{
            *f=0;
        }
        return 0;
}

// integrations
double ComputeSigmaR::LNLOdip(double Q, double x) { // old LNLObeufDIP
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOdip_z2(double Q, double x) { // old LNLObeufDIP
    double integral, error, prob;
    const int ndim=3;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip_z2,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsub(double Q, double x) { // old LNLObeufQGiancu
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgunsub,&userdata,&integral,&error,&prob);
    // return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
    return 2*fac*2.0*M_PI*integral; // removed x01, x02 integral jacobians --- use with logarithmic x01, x02 integrals
}

double ComputeSigmaR::LNLOsigma3(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILsigma3,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgsub(double Q, double x) { // no old implementation/version. CXY is similar but distinct, sub vs. xbj-sub.
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgsub,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgunsubRisto,&userdata,&integral,&error,&prob);
    return (2/2)*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral; // TODO a difference of a factor of 2 w.r.t. BEUF formulation was discovered in initial testing
}

double ComputeSigmaR::LNLOqgsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgsubRisto,&userdata,&integral,&error,&prob);
    return (2/2)*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral; // TODO a difference of a factor of 2 w.r.t. BEUF formulation was discovered in initial testing
}
//*/



///===========================================================================================
/// --- LLL ---  NLO  MASSIVE --- LLL ---------------------

int integrand_ILdip_massive_LiLogConst(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double res;

    res = SKernel*(ILdip_massive_LiLogConst(Q,z1,x01sq,mf))*x01*alphfac;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILdip_massive_Iab(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);
    double xi=x[2]; // New integration variable in I_ab terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double res;

    res = SKernel*(ILdip_massive_Iab(Q,z1,x01sq,mf,xi))*x01*alphfac;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILdip_massive_Icd(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);
    double xi=x[2]; // 1st new integration variable in I_cd terms
    double intx=x[3]; // 2nd new integration variable in I_cd terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double res;

    res = SKernel*(ILdip_massive_Icd(Q,z1,x01sq,mf,xi,intx))*x01*alphfac;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsub_massive(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_u1 = x[6];
    double y_t2 = x[7];
    double y_u2 = x[8];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    // if (isnan(x01) or isnan(x02) or isnan(sqrt(x21sq))){
    //     cout << x01 << " " << Xrpdt << " " << x02 << " " << Xrpdt << " " << sqrt(x21sq) << " " << Xrpdt << endl; 
    //     cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << x[8] << endl;
    //     *f=0;
    //     return 0;
    // }
    double SKernel_dipole = 1.0 - Optr->Sr(x01,Xrpdt);
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);


    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double dipole_term  = SKernel_dipole  * ILNLOqg_massive_dipole_part_unintegrated(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_u1,y_t2,y_u2); // Terms proportional to N_01
    double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_unintegrated(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_u1,y_t2,y_u2); // Terms proportional to N_012

    double res =   jac*alphfac*( dipole_term + tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsub_massive_I1(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;
    double SKernel_dipole = 1.0 - Optr->Sr(x01,Xrpdt);
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double dipole_term  = SKernel_dipole  * ILNLOqg_massive_dipole_part_I1(Q,mf,z1,z2,x01sq,x02sq,x21sq); // Terms proportional to N_01
    double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I1(Q,mf,z1,z2,x01sq,x02sq,x21sq); // Terms proportional to N_012

    double res =   jac*alphfac*( dipole_term + tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsub_massive_I2(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_u1 = x[6]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I2(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_u1); // Terms proportional to N_012
    // double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I2_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsub_massive_I2_fast(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    // double y_u1 = x[6]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I2_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsub_massive_I3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_t2 = x[6];
    double y_u1 = x[7]; // not used in fast
    double y_u2 = x[8]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    // if (isnan(x01) or isnan(x02) or isnan(sqrt(x21sq))){
    //     cout << x01 << " " << Xrpdt << " " << x02 << " " << Xrpdt << " " << sqrt(x21sq) << " " << Xrpdt << endl; 
    //     cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << x[8] << endl;
    //     *f=0;
    //     return 0;
    // }
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);


    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I3(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_u1,y_t2,y_u2); // Terms proportional to N_012
    // double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I3_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_t2); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ILqgunsub_massive_I3_fast(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_t2 = x[6];
    // double y_u1 = x[7]; // not used in fast
    // double y_u2 = x[8]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);


    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ILNLOqg_massive_tripole_part_I3_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_t2); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}



double ComputeSigmaR::LNLOdip_massive_LiLogConst(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip_massive_LiLogConst,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}


double ComputeSigmaR::LNLOdip_massive_Iab(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=3; // One more integral than in the Li & Log & Const terms
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip_massive_Iab,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOdip_massive_Icd(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=4; // Two more integrals than in the Li & Log & Const terms
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILdip_massive_Icd,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsub_massive(double Q, double x, double mf) {
    double integral, error, prob;
    // const int ndim=5; // 5 + 2*2 deterministic cuhre dimensions from G_qg*G_qg evaluated inside the integral
    const int ndim=9; // MC full phase space
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgunsub_massive,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsub_massive_I1(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=5; 
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ILqgunsub_massive_I1,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsub_massive_I2(double Q, double x, double mf) {
    double integral, error, prob;
    int ndim;
    // const int ndim=7; // MC 5+2 phase space
    // const int ndim=6; // fast MC 5+1 phase space
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    if (nlodis_config::PERF_MODE == nlodis_config::DISABLED) {
        ndim = 7;
        Cuba(cubamethod,ndim,integrand_ILqgunsub_massive_I2,&userdata,&integral,&error,&prob);
    } else if (nlodis_config::PERF_MODE == nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION) {
        ndim = 6;
        Cuba(cubamethod,ndim,integrand_ILqgunsub_massive_I2_fast,&userdata,&integral,&error,&prob);
    }
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::LNLOqgunsub_massive_I3(double Q, double x, double mf) {
    double integral, error, prob;
    int ndim;
    // const int ndim=9; // MC full 5+2+2 phase space
    // const int ndim=7; // fast MC full 5+2 phase space
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    if (nlodis_config::PERF_MODE == nlodis_config::DISABLED) {
        ndim = 9;
        Cuba(cubamethod,ndim,integrand_ILqgunsub_massive_I3,&userdata,&integral,&error,&prob);
    } else if (nlodis_config::PERF_MODE == nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION) {
        ndim = 7;
        Cuba(cubamethod,ndim,integrand_ILqgunsub_massive_I3_fast,&userdata,&integral,&error,&prob);
    }
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}



///===========================================================================================
///===========================================================================================
///===========================================================================================
/*
// --- T T T --- NLO -------- T T T --- NLO -------- T T T --- NLO -----------
*/

///*
// NLO
double ComputeSigmaR::ITdip(double Q, double z1, double x01sq) {
    double fac1 = z1*(1-z1)*(Sq(z1)+Sq(1.0 - z1));
    double facNLO = Sq(Q*gsl_sf_bessel_K1(Q*sqrt(z1*(1-z1)*x01sq)));
    return fac1*facNLO;
}

int integrand_ITdip(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) { // integrand_ITbeufDIP
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double regconst = 5.0/2.0 - Sq(M_PI)/6.0;
    double res;

    res = SKernel*( Optr->ITdip(Q,z1,x01sq) )*x01*( alphfac*(0.5*Sq(log(z1/(1-z1))) + regconst ) );

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITdip_z2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    double res;
    if (z2min > 0.5){ // Check that z2min is not too large since z1 limits are [z2min , 1-z2min]. IF it is too large, return res=0.
        res = 0;
        *f=res;
    return 0;
    }
    double z1=((1.0-z2min)-z2min)*x[0]+z2min;
    double z2=((z1)-z2min)*x[1]+z2min; // z2 integration variable for k1-term AND for the mirror symmetry combined integral!
    double x01=nlodis_config::MAXR*x[2];
    double x01sq=Sq(x01);
    double jac=((1.0-z2min)-z2min)*((z1)-z2min); // jacobiaaniin skaalaus molemmista z-muuttujista!
    double Xrpdt= z2min * X0/z2; // consistent with the qg-terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdt); //ok?

    double fivefourths = 5.0/4.0;
    double loopcontribution = 1.0/z2 * ( log(1.0 + z2 / (1.0-z1) ) + log(1.0 - z2 / z1) ) + 1.0/z1 * fivefourths ; // a.k.a k1- or z1-term
    res = 2 * jac * SKernel*(Optr->ITdip(Q,z1,x01sq))*x01*( alphfac*( loopcontribution ) ); // factor of 2 is from the combination of the mirror symmteric q/qbar terms.

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double x01=sqrt(x01sq);
    double x02=sqrt(x02sq);
    double x21=sqrt(x21sq);

    double X3sq = z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;

    double Qsq = Sq(Q);
    double x_x01 = cov_to_eta_x(x01, x, Qsq);
    double x_x02 = cov_to_eta_x(x02, x, Qsq);
    double x_x21 = cov_to_eta_x(x21, x, Qsq);

    double bessel_innerfun = Q*sqrt(X3sq);
    double facNLO = 0;
    if (bessel_innerfun > 200){
        // cout << "bessel_innerfun = " << bessel_innerfun << " Q " << Q << " X3sq " << X3sq << endl;
        return 0;
    } else if (bessel_innerfun < 1e-30){
        return 1;
    }else{
        facNLO = Sq(Q*gsl_sf_bessel_K1(bessel_innerfun));
    }
    double res = facNLO*(1-SrTripole(x01, x_x01, x02, x_x02, x21, x_x21));
    return res;
}

double ComputeSigmaR::ITNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) { // old ITNLObeufQG
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double X3sq 	= z1*(1.0 - z1 - z2)*x01sq + z2*(1.0 - z1 - z2)*x02sq + z2*z1*x21sq;

    double fac1 	= z1*(1.0 - z1);
    double xi 		= z2/(1.0-z1);
    double fac11	= Sq(z1)+Sq(1-z1);
    double fac12 	= 1+Sq(1-xi);
    double fun2 	= x20x21/(x02sq*x21sq);
    double fac13 	= 1/x02sq - fun2;
    double fac3 	= Sq(xi)*(	fac11*fun2	+ 2*fac1*(1-xi)*x20x21/(x02sq*X3sq)	-	(1-z1)*(1-xi)*(z1+xi-z1*xi)/X3sq	);

    double facNLO1 = Bessel1Tripole(Q, x, z1, z2, x01sq, x02sq, x21sq);
    double facNLO2 = Bessel1Tripole(Q, x, z1, z2, x01sq, 0    , x01sq);

    double res = fac1*((fac11*fac12*fac13)*(facNLO1 - facNLO2) + fac3*facNLO1);
    return res;
}

double ComputeSigmaR::ITNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) {
    double impactfac_lo = ITLO(Q,z1,x01sq);
    double BKkernel = (K_kernel(x01sq,x02sq,x21sq) - 1.0)*2*K_lobk(x01sq,x02sq,x21sq);
    double S012 = SrTripole(sqrt(x01sq),x,sqrt(x02sq),x,sqrt(x21sq),x);
    double Sdipole_bkevol = Sr(sqrt(x01sq),x) - S012;

    double res = BKkernel*Sdipole_bkevol*impactfac_lo;
    return res;
}

double ComputeSigmaR::ITNLOqgRisto(double Q, double x, double z0, double z2, double x01sq, double x02sq, double x21sq){
    double x20x21 = -0.5*(x01sq - x21sq - x02sq);
    double z3 = 1 - z0 - z2;
    double X3sq = z0*z3*x01sq + z0*z2*x02sq + z2*z3*x21sq;

    double fun2 = x20x21/(x02sq*x21sq);

    double factermF = z3/(z0+z2) * ( 1-2*z3*(1-z3) ) * ( 2*z0*(z0+z2) + Sq(z2) )/x02sq;
    double factermG = z0/(z2+z3) * ( 1-2*z0*(1-z0) ) * ( 2*z3*(z2+z3) + Sq(z2) )/x21sq;
    double factermPI = (-2*z0*z3*(
                        (Sq(1-z0)+Sq(z0)) + (Sq(1-z3)+Sq(z3))
                        )*X3sq*fun2
                        +(2*z0*Sq(z2*z3))/(z0+z2)*x20x21/(x02sq)
                        +(2*z3*Sq(z2*z0))/(z2+z3)*x20x21/(x21sq)
                        +z0*z2*z3*(
                            Sq(z0+z2) + Sq(z2+z3) +2*(Sq(z0) + Sq(z3))
                            +Sq(z0*z3)/(Sq(z0+z2))+Sq(z0*z3)/(Sq(z2+z3))
                            -z2*( (z0+z2)/(z2+z3)+ (z2+z3)/(z0+z2) )
                            -( (1-2*z3*(1-z3))*(2*z0*(z0+z2)+Sq(z2)) )/(Sq(z0+z2))
                            -( (1-2*z0*(1-z0))*(2*z3*(z2+z3)+Sq(z2)) )/(Sq(z2+z3))
                            )
                        )/X3sq;

    double facBesselTrip = Bessel1Tripole(Q, x, z3, z2, x01sq, x02sq, x21sq);

    double Qbarn = sqrt(z3*(1-z3))*Q;
    double Qbaro = sqrt(z0*(1-z0))*Q;
    double SKernel = 1.0 - Sr(sqrt(x01sq),x);

    ///*
    double facExpn, facExpo;
    if(x02sq/x01sq>1e-8 && x02sq/x01sq<5e2){
        facExpn = std::exp(-x02sq/(std::exp(M_EULER)*x01sq));
    }else if (x02sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }
    if(x21sq/x01sq>1e-8 && x21sq/x01sq<5e2){
        facExpo = std::exp(-x21sq/(std::exp(M_EULER)*x01sq));
    }else if (x21sq/x01sq<1e-8){
        facExpn = 1;
    }else{ facExpn = 0; }

    double facBesExpn, facBesExpo;
    if(Qbarn*sqrt(x01sq)>1e-8 && Qbarn*sqrt(x01sq)<5e2){
        facBesExpn = Sq(Q*gsl_sf_bessel_K1(Qbarn*sqrt(x01sq)))*facExpn*SKernel;
    }else{ facBesExpn = 0; }
    if(Qbaro*sqrt(x01sq)>1e-8 && Qbaro*sqrt(x01sq)<5e2){
        facBesExpo = Sq(Q*gsl_sf_bessel_K1(Qbaro*sqrt(x01sq)))*facExpo*SKernel;
    }else{ facBesExpo = 0; }
    //*/
    //double facBesExpn = Sq(Q*gsl_sf_bessel_K1(Qbarn*sqrt(x01sq)))*std::exp(-x02sq/(std::exp(M_EULER)*x01sq))*SKernel;
    //double facBesExpo = Sq(Q*gsl_sf_bessel_K1(Qbaro*sqrt(x01sq)))*std::exp(-x21sq/(std::exp(M_EULER)*x01sq))*SKernel;

    double term1 = factermF*( facBesselTrip - facBesExpn );
    double term2 = factermG*( facBesselTrip - facBesExpo );
    double term3 = factermPI*facBesselTrip;
    double res = term1 + term2 + term3;

    return res;
}

int integrand_ITqgunsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) { // old integrand_ITbeufQGian
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min; // linear z_2 int
    // double ymax = std::log(1-z1);
    // double ymin = std::log(z2min);
    // double y=(ymax-ymin)*x[1]+ymin;
    // double z2=std::exp(y);
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    // double rhomax=std::log(nlodis_config::MAXR);
    // double rhomin=std::log(nlodis_config::MINR);
    // double rho01=(rhomax-rhomin)*x[2] + rhomin;
    // double rho02=(rhomax-rhomin)*x[3] + rhomin;
    // double x01=std::exp(rho01);
    // double x02=std::exp(rho02);
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min)*Sq(nlodis_config::MAXR); // linear z_2 int, linear in x01 x02
    // double jac=(1.0-z2min)*(ymax-ymin)*z2; // log z_2 int
    // double jac=(1.0-z2min)*(1.0-z1-z2min)*Sq(rhomax-rhomin)*x01*x02; // linear z_2 int, log x01, log x02
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res = jac*alphfac*( Optr->ITNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }
    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double res =   jac*alphfac*( Optr->ITNLOsigma3(Q,Xrpdt,z1,0,x01sq,x02sq,x21sq) )/z2*x01*x02; // K_NLO(z_2 = 0) kirjoitettu explisiittisesti

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min >= 1.0){ // Check that z2min is not too large. IF it is too large, return res=0.
        *f=0;
        return 0;
    }
    double z1=x[0];
    double z2=(1.0-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double ITNLOqg_z2 = 0;
    if (1-z1-z2 > 0){ // heaviside theta front factor on the full NLO term
        ITNLOqg_z2 = (Optr->ITNLOqg(Q,Xrpdt,z1,z2,x01sq,x02sq,x21sq));
    }
    double ITNLOqg_z2_to_0 = Optr->ITNLOqg_subterm(Q,Xrpdt,z1,0,x01sq,x02sq,x21sq);

    double res = jac*alphfac*(
                    ITNLOqg_z2
                    -
                    ITNLOqg_z2_to_0
                )/(z2)*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=(1.0-x0lim)*x[0];
    double z2=(1.0-z0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim)*(1.0-z0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;

    double res;
    res = jac*alphfac*( Optr->ITNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq) )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata){
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0; //0.01;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double x0lim=xbj/X0;
    double z0=x[0];
    double z2=(1.0-x0lim)*x[1]+x0lim;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-x0lim);
    double Xrpdt=xbj/z2;

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;

    double res;
    res =   jac*alphfac*(
                    (Optr->ITNLOqgRisto(Q,Xrpdt,z0,z2,x01sq,x02sq,x21sq))*(Optr->heaviside_theta(1-z0-z2))
                    -Optr->ITNLOqgRisto(Q,Xrpdt,z0,0,x01sq,x02sq,x21sq)
                )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}


// integrations
double ComputeSigmaR::TNLOdip(double Q, double x) { // old name: TNLObeufDIP
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOdip_z2(double Q, double x) { // old name: TNLObeufDIP
    double integral, error, prob;
    const int ndim=3;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip_z2,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgunsub(double Q, double x) { // old TNLObeufQGiancu
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgunsub,&userdata,&integral,&error,&prob);
    // return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
    return 2*fac*2.0*M_PI*integral; // removed x01 x02 integral jacobian --- use with logarithmic x01 x02 integral
}

double ComputeSigmaR::TNLOsigma3(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITsigma3,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgsub(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgsub,&userdata,&integral,&error,&prob);
    return 2*fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgunsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgunsubRisto,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgsubRisto(double Q, double x) {
    double integral, error, prob;
    const int ndim=5;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef;
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgsubRisto,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}






///===========================================================================================
/// --- TTT ---  NLO  MASSIVE --- TTT ---------------------


int integrand_ITdip_massive_I1(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double res;

    res = SKernel*(ITdip_massive_0(Q,z1,x01sq,mf))*x01*alphfac;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITdip_massive_I2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);
    double xi=x[2]; // New integration variable in I_ab terms

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double res;

    res = SKernel*(ITdip_massive_1(Q,z1,x01sq,mf,xi))*x01*alphfac;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITdip_massive_I3(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z1=x[0];
    double x01=nlodis_config::MAXR*x[1];
    double x01sq=Sq(x01);
    double y_chi=x[2];
    double y_u=x[3];

    double alphabar=Optr->Alphabar(x01sq);
    double alphfac=alphabar*CF/Nc;
    double Xrpdty_lo = Optr->Xrpdty_DIP(xbj, Sq(Q), x01sq);
    double SKernel = 1.0 - Optr->Sr(x01,Xrpdty_lo);
    double res;

    res = SKernel*(ITdip_massive_2(Q,z1,x01sq,mf,y_chi,y_u))*x01*alphfac;
    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsub_massive_I1(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;
    double SKernel_dipole = 1.0 - Optr->Sr(x01,Xrpdt);
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double dipole_term  = SKernel_dipole  * ITNLOqg_massive_dipole_part_I1(Q,mf,z1,z2,x01sq,x02sq,x21sq); // Terms proportional to N_01
    double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I1(Q,mf,z1,z2,x01sq,x02sq,x21sq); // Terms proportional to N_012

    double res =   jac*alphfac*( dipole_term + tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsub_massive_I2(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_u1 = x[6]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I2(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_u1); // Terms proportional to N_012
    // double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I2_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsub_massive_I2_fast(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    // double y_u1 = x[6]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;
    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);

    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I2_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsub_massive_I3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_t2 = x[6];
    double y_u1 = x[7]; // not used in fast
    double y_u2 = x[8]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=std::max(0.0, x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102));
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);


    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I3(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_u1,y_t2,y_u2); // Terms proportional to N_012
    // double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I3_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_t2); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

int integrand_ITqgunsub_massive_I3_fast(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) {
    Userdata *dataptr = (Userdata*)userdata;
    double Q=dataptr->Q;
    double xbj=dataptr->xbj;
    double X0=dataptr->icX0;
    double mf=dataptr->qMass;
    ComputeSigmaR *Optr = dataptr->ComputerPtr;
    double z2min = Optr->z2lower_bound(xbj,Sq(Q));
    if (z2min > 1.0){ // Check that z2min is not too large. IF it is too large, return *f=0.
        *f=0;
        return 0;
    }

    double z1=(1.0-z2min)*x[0];
    double z2=((1.0-z1)-z2min)*x[1]+z2min;
    double x01=nlodis_config::MAXR*x[2];
    double x02=nlodis_config::MAXR*x[3];
    double phix0102=2.0*M_PI*x[4];
    double y_t1 = x[5];
    double y_t2 = x[6];
    // double y_u1 = x[7]; // not used in fast
    // double y_u2 = x[8]; // not used in fast
    double x01sq=Sq(x01);
    double x02sq=Sq(x02);
    double x21sq=x01sq+x02sq-2.0*sqrt(x01sq*x02sq)*cos(phix0102);
    double jac=(1.0-z2min)*(1.0-z1-z2min);
    double Xrpdt= Optr->Xrpdty_NLO(Q*Q, z2, z2min, X0, x01sq, x02sq, x21sq); //z2min * X0/z2;

    double SKernel_tripole = 1.0 - Optr->SrTripole(x01,Xrpdt,x02,Xrpdt,sqrt(x21sq),Xrpdt);


    Alphasdata alphasdata;
    alphasdata.x01sq=x01sq;
    alphasdata.x02sq=x02sq;
    alphasdata.x21sq=x21sq;
    double alphabar=Optr->Alphabar_QG( &alphasdata );
    double alphfac=alphabar*CF/Nc;

    double tripole_term = SKernel_tripole * ITNLOqg_massive_tripole_part_I3_fast(Q,mf,z1,z2,x01sq,x02sq,x21sq,y_t1,y_t2); // Terms proportional to N_012

    double res =   jac*alphfac*( tripole_term )/z2*x01*x02;

    if(gsl_finite(res)==1){
        *f=res;
    }else{
        *f=0;
    }
    return 0;
}

double ComputeSigmaR::TNLOdip_massive_I1(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=2;
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip_massive_I1,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}


double ComputeSigmaR::TNLOdip_massive_I2(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=3; // One more integral than in the Li & Log & Const terms
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip_massive_I2,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOdip_massive_I3(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=4; // Two more integrals than in the Li & Log & Const terms
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITdip_massive_I3,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*integral;
}


// ----- qqbarg contribution -----
double ComputeSigmaR::TNLOqgunsub_massive_I1(double Q, double x, double mf) {
    double integral, error, prob;
    const int ndim=5; 
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    Cuba(cubamethod,ndim,integrand_ITqgunsub_massive_I1,&userdata,&integral,&error,&prob);
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgunsub_massive_I2(double Q, double x, double mf) {
    double integral, error, prob;
    int ndim;
    // const int ndim=7; // MC 5+2 phase space
    // const int ndim=6; // fast MC 5+1 phase space
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    if (nlodis_config::PERF_MODE == nlodis_config::DISABLED) {
        ndim = 7;
        Cuba(cubamethod,ndim,integrand_ITqgunsub_massive_I2,&userdata,&integral,&error,&prob);
    } else if (nlodis_config::PERF_MODE == nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION) {
        ndim = 6;
        Cuba(cubamethod,ndim,integrand_ITqgunsub_massive_I2_fast,&userdata,&integral,&error,&prob);
    }
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}

double ComputeSigmaR::TNLOqgunsub_massive_I3(double Q, double x, double mf) {
    double integral, error, prob;
    int ndim;
    // const int ndim=9; // MC full 5+2+2 phase space
    // const int ndim=7; // fast MC full 5+2 phase space
    double fac=4.0*Nc*alphaem/Sq(2.0*M_PI)*sumef_from_mass(mf);
    Userdata userdata;
    userdata.Q=Q;
    userdata.xbj=x;
    userdata.icX0=icX0;
    userdata.qMass=mf;
    userdata.ComputerPtr=this;
    if (nlodis_config::PERF_MODE == nlodis_config::DISABLED) {
        ndim = 9;
        Cuba(cubamethod,ndim,integrand_ITqgunsub_massive_I3,&userdata,&integral,&error,&prob);
    } else if (nlodis_config::PERF_MODE == nlodis_config::MASSIVE_EXPLICIT_BESSEL_DIM_REDUCTION) {
        ndim = 7;
        Cuba(cubamethod,ndim,integrand_ITqgunsub_massive_I3_fast,&userdata,&integral,&error,&prob);
    }
    return fac*2.0*M_PI*nlodis_config::MAXR*nlodis_config::MAXR*integral;
}


//*/
