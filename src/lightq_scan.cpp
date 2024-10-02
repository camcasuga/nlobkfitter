#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <csignal>
#include <ctime>
#include <gsl/gsl_errno.h>
#include <tools/tools.hpp>
#include <sstream>
#include <unistd.h>

#include <amplitudelib/amplitudelib.hpp>
#include <tools/interpolation.hpp>
#include "solver.hpp"
#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"

#include "data.hpp"
#include "nlodissigmar.hpp"
#include "nlodis_config.hpp"

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnScan.h>
#include <Minuit2/MnMinimize.h>
#include <Minuit2/MnSimplex.h>

 #include <gsl/gsl_errno.h>

using namespace std;
using namespace ROOT::Minuit2;
//void ErrHandler(const char * reason,const char * file,int line,int gsl_errno);
int errors_mmyiss;
void ErrHandlerCustom(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno)
{   // 14 = failed to reach tolerance
    // 18 = roundoff error prevents tolerance from being achieved
    // 11 = maximum number of subdivisions reached
    // 15: underflows

    if (gsl_errno == 11) return; // ignore max subdivision errors
    if (gsl_errno == 18) return; // roundoff errors from small r?
    // if (gsl_errno == 15 or gsl_errno == 16) return;
    // Ugly hack, comes from the edges of the z integral in virtual_photon.cpp
    // Overflows come from IPsat::bint when it is done analytically
    // Hope is that these errors are handled correctly everywhere
    errors_mmyiss++;
    std::cerr << file << ":"<< line <<": Error " << errors_mmyiss << ": " <<reason
            << " (code " << gsl_errno << ")." << std::endl;
}

int main( int argc, char* argv[] )
{
    gsl_set_error_handler(&ErrHandlerCustom);

    // NLO DIS SIGMA_R COMPUTATION CONFIGS
//     nlodis_config::CUBA_EPSREL = 10e-3;
    nlodis_config::CUBA_EPSREL = 2e-3; // highacc def1
//     nlodis_config::CUBA_MAXEVAL= 2e7;
    nlodis_config::CUBA_MAXEVAL= 5e7; // highacc def1
    nlodis_config::MINR = 1e-6;
    nlodis_config::MAXR = 30;
    nlodis_config::PRINTDATA = true;
    bool useNLO = true;
    bool computeNLO = useNLO;
    nlodis_config::USE_MASSES = false;

    string cubaMethod = "vegas";
    // string cubaMethod = "suave";

    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::KINEMATICAL_CONSTRAINT = config::KC_NONE;

    config::VERBOSE = true;
    //config::RINTPOINTS = 512/4;
    //config::THETAINTPOINTS = 512/4;

    // config::INTACCURACY = 10e-3;//0.02;
    config::INTACCURACY = 5e-3; // highacc def1
    // config::INTACCURACY = 20e-3; // quick low acc
    config::MINR = 1e-6;
    config::MAXR = 30;
    config::RPOINTS = 100;
    config::DE_SOLVER_STEP = 0.4; // Rungekutta step

    // Constants
    config::NF=3;   // Only light quarks
    config::LAMBDAQCD = 0.241;


    Data data;
    data.SetMinQsqr(0.75);
    data.SetMaxQsqr(50);
    data.SetMaxX(0.01);
    data.LoadData("./data/light_quark_f2/hera_I_combined_eplus_lightq", TOTAL);

    MnUserParameters parameters;

    bool useSUB, useResumBK, useKCBK, useImprovedZ2Bound, useBoundLoop;
    bool useSigma3 = false;
    string helpstring = "Argument order: SCHEME BK RC useImprovedZ2Bound useBoundLoop [Qs0 C^2 gamma] X0_if X0_bk e_c Q0sq Y0 eta0\nsub/unsub/unsub+ resumbk/trbk/lobk parentrc/guillaumerc/fixedrc z2improved/z2simple z2boundloop/unboundloop";
    string string_sub, string_bk, string_rc;
    if (argc<2){ cout << helpstring << endl; return 0;}
    // Argv[0] is the name of the program

    string_sub = string(argv [1]);
    if (string(argv [1]) == "sub"){
        useSUB = true;
        nlodis_config::SUB_SCHEME = nlodis_config::SUBTRACTED;
    } else if (string(argv [1]) == "unsub"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = false;
    } else if (string(argv [1]) == "unsub+"){
        nlodis_config::SUB_SCHEME = nlodis_config::UNSUBTRACTED;
        useSUB = false;
        useSigma3 = true;
    } else {cout << helpstring << endl; return -1;}

    string_bk = string(argv [2]);
    if (string(argv [2]) == "resumbk"){
            config::EULER_METHOD = false;    // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = true;       // Resum doulbe logs
            config::RESUM_SINGLE_LOG = true; // Resum single logs
            config::KSUB = 0.65;             // Optimal value for K_sub
            config::NO_K2 = true;            // Do not include numerically demanding full NLO part
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_RESUM;
    }else if (string(argv [2]) == "kcbk"){
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_BEUF_K_PLUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_KCBK_BEUF;
    }else if (string(argv [2]) == "trbk"){  // Target Rapidity BK
            config::EULER_METHOD = true;     // Kinematical constraint requires this
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_EDMOND_K_MINUS;
            config::DE_SOLVER_STEP = 0.05;  //0.02; // Euler method requires smaller step than RungeKutta!
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_TRBK_EDMOND;
            nlodis_config::TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_RQ0;
            // nlodis_config::TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_QQ0;
    }else if (string(argv [2]) == "lobk"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_EXPLICIT;
    }else if (string(argv [2]) == "lobkold"){
            config::EULER_METHOD = false;   // Use Runge-Kutta since no kin. constraint
            config::RESUM_DLOG = false;
            config::RESUM_SINGLE_LOG = false;
            config::KINEMATICAL_CONSTRAINT = config::KC_NONE;
            nlodis_config::SUB_TERM_KERNEL = nlodis_config::SUBTERM_LOBK_Z2TOZERO;
    } else {cout << helpstring << endl; return -1;}

    string_rc = string(argv [3]);
    if (string(argv [3]) == "parentrc" or string(argv [3]) == "pdrc"){
            config::RC_LO = config::PARENT_LO;
            config::RESUM_RC = config::RESUM_RC_PARENT;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_PARENT;
    } else if (string(argv [3]) == "guillaumerc" or string(argv [3]) == "gbrc"){
            config::RC_LO = config::GUILLAUME_LO;
            config::RESUM_RC = config::RESUM_RC_GUILLAUME;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_GUILLAUME;
    } else if (string(argv [3]) == "fixedrc" or string(argv [3]) == "fc"){
            config::RC_LO = config::FIXED_LO;
            config::RESUM_RC = config::RESUM_RC_FIXED;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_FIXED;
    } else if (string(argv [3]) == "smallestrc" or string(argv [3]) == "sdrc"){
            config::RC_LO = config::SMALLEST_LO;
            config::RESUM_RC = config::RESUM_RC_SMALLEST;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_SMALLEST;
    } else if (string(argv [3]) == "balitskysmallrc" or string(argv [3]) == "balsdrc"){
            // Even though this coupling should be more realistic than smallest dipole alone,
            // the subtraction between the LO and qg terms is not exact. This shortcoming makes
            // this coupling less than ideal.
            config::RC_LO = config::BALITSKY_LO;
            config::RESUM_RC = config::RESUM_RC_SMALLEST;
            nlodis_config::RC_DIS = nlodis_config::DIS_RC_SMALLEST;
    }else {cout << helpstring << endl; return -1;}

    if (string(argv [4]) == "z2improved" or string(argv [4]) == "z2imp"){
        nlodis_config::Z2MINIMUM = nlodis_config::Z2IMPROVED;
        useImprovedZ2Bound = true;
    } else if (string(argv [4]) == "z2simple" or string(argv [4]) == "z2sim"){
        nlodis_config::Z2MINIMUM = nlodis_config::Z2SIMPLE;
        useImprovedZ2Bound = false;
    } else {cout << helpstring << endl; return -1;}

    if (string(argv [5]) == "z2boundloop" or string(argv [5]) == "z2b"){
        useBoundLoop = true;
    } else if (string(argv [5]) == "unboundloop" or string(argv [5]) == "unb"){
        useBoundLoop = false;
    } else {cout << helpstring << endl; return -1;}

    int argi=6;
    double icqs0sq, iccsq, icx0_if, icx0_bk, ic_ec, icgamma, icQ0sq, icY0, icEta0;
    // reading fit initial condition parameters:
    if (argc == argi){
        cout << "Swarmscan needs Q^2, C^2 and gamma IC values at the least!" << endl;
        cout << helpstring << endl;
        exit(1);
    } else if (argc == argi+3){
        icqs0sq   = stod( argv [argi] ); argi++;
        iccsq     = stod( argv [argi] ); argi++;
        icgamma   = stod( argv [argi] ); argi++;
        icx0_if   = 1.0;
        icx0_bk   = 0.01;
        ic_ec     = 1.0;
        icQ0sq    = 1.0;
        icY0      = 0;
        icEta0    = 0;
    } else {
        icqs0sq   = stod( argv [argi] ); argi++;
        iccsq     = stod( argv [argi] ); argi++;
        icgamma   = stod( argv [argi] ); argi++;
        icx0_if   = stod( argv [argi] ); argi++;
        icx0_bk   = stod( argv [argi] ); argi++;
        ic_ec     = stod( argv [argi] ); argi++;
        icQ0sq    = stod( argv [argi] ); argi++;
        icY0      = stod( argv [argi] ); argi++;
        icEta0    = stod( argv [argi] ); argi++;
    }


    // Constructing Q and C^2 grids
    double Qcalc = icqs0sq/1000.; // input in milli-GeV² integer strings

    double Cdown = 0.1;
    double Cup = 5.0;
    double CStep = pow(Cup/Cdown , iccsq/240.);
    double Ccalc = Cdown * CStep; // compute C² on a logarithmic grid

    double old_sigma02 = 1.;
    double mass_charm = 1.27; // MSbar value as default
    cout << Qcalc << " " << Ccalc << endl;

	parameters.Add("qs0sqr",		        Qcalc);
    // parameters.Add("fitsigma0",		        1.0); // 1mb = 2.568 GeV²
    parameters.Add("alphascalingC2",	    Ccalc);
    parameters.Add("p",                   ic_ec );
    //parameters.Add("p",                     ic_p);
    parameters.Add("anomalous_dimension",   icgamma );
    parameters.Add("icx0_nlo_impfac",       icx0_if );
    parameters.Add("icx0_bk",               icx0_bk );
    parameters.Add("icTypicalPartonVirtualityQ0sqr", icQ0sq );
    parameters.Add("initialconditionY0",    icY0 );
    parameters.Add("eta0", icEta0 );
    parameters.Add("old_sigma02", old_sigma02 );
    parameters.Add("mass_charm", mass_charm );

    NLODISFitter fitter(parameters);
    fitter.AddDataset(data);
    fitter.SetNLO(useNLO);
    fitter.SetSUB(useSUB);
    fitter.SetSigma3(useSigma3);
    fitter.UseImprovedZ2Bound(useImprovedZ2Bound);
    fitter.UseConsistentlyBoundLoopTerm(useBoundLoop);
    fitter.SetCubaMethod(cubaMethod);

    cout << "=== Perturbative settings ===" << endl;
    cout << std::boolalpha;
    cout    << "# === Perturbative settings ===" << endl
            << "# Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "# Use LOBK (DL,SL==false): " << (!(config::RESUM_DLOG) 
                                    and !(config::RESUM_SINGLE_LOG)) << endl
            << "# Use ResumBK (DL,SL==true,KC_NONE): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == true)) << endl
            << "# Use NLOBK (DL,SL==true,KC_NONE,NO_K2==false): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == false)) << endl
            << "# KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << " (0 BEUF_K_PLUS, 1 EDMOND_K_MINUS, 2 NONE)" << endl
            << "# Target eta0 RHO shift: " << nlodis_config::TRBK_RHO_PRESC << " (0 TRBK_RHO_DISABLED, 1 TRBK_RHO_QQ0, 2 TRBK_RHO_RQ0)" << endl
            << "# Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 2 parent_beta, 3 smallest, 4 balitsky, 5 frac, 6 guillaume)" << endl
            << "# Running Coupling: (RC_NLO):    " << config::RC_NLO << " (0 fc, 1 parent, 2 smallest)" << endl
            << "# Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 fc, 1 balitsky, 2 parent, 3 smallest, 4 guillaume)" << endl
            << "# Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 smallest, 3 guillaume)" << endl
            << "# Use NLOimpact: " << useNLO << endl
            << "# Use SUBscheme: " << useSUB << endl
            << "# Use Sigma3: " << useSigma3 << endl
            << "# Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "# Use Z2 loop term: " << useBoundLoop << endl
            << "# Cuba MC: " << cubaMethod
                << ", Cuba eps = " << nlodis_config::CUBA_EPSREL
                << ", Cuba maxeval = " << (float)nlodis_config::CUBA_MAXEVAL
                << endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;
    cout << "=== Initial parameters ===" << endl;
    cout << parameters << endl;
    if(nlodis_config::VERBOSE) cout << "=== Starting fit ===" << endl;

    MnMinimize fit(fitter, parameters);
    FunctionMinimum min = fit();
    std::cout<<"minimum: "<<min<<std::endl;

    cout << "=== Perturbative settings were ===" << endl;
    cout << std::boolalpha;
    cout    << "# === Perturbative settings ===" << endl
            << "# Settings: " << string_sub << " (scheme), " << string_bk << ", " << string_rc << endl
            << "# Use LOBK (DL,SL==false): " << (!(config::RESUM_DLOG) 
                                    and !(config::RESUM_SINGLE_LOG)) << endl
            << "# Use ResumBK (DL,SL==true,KC_NONE): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == true)) << endl
            << "# Use NLOBK (DL,SL==true,KC_NONE,NO_K2==false): " << ((config::RESUM_DLOG) 
                                    and (config::RESUM_SINGLE_LOG)
                                    and (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
                                    and (config::NO_K2 == false)) << endl
            << "# KinematicalConstraint / target eta0 BK: " << config::KINEMATICAL_CONSTRAINT << " (0 BEUF_K_PLUS, 1 EDMOND_K_MINUS, 2 NONE)" << endl
            << "# Target eta0 RHO shift: " << nlodis_config::TRBK_RHO_PRESC << " (0 TRBK_RHO_DISABLED, 1 TRBK_RHO_QQ0, 2 TRBK_RHO_RQ0)" << endl
            << "# Running Coupling: (RC_LO):    " << config::RC_LO << " (0 fc, 1 parent, 2 parent_beta, 3 smallest, 4 balitsky, 5 frac, 6 guillaume)" << endl
            << "# Running Coupling: (RC_NLO):    " << config::RC_NLO << " (0 fc, 1 parent, 2 smallest)" << endl
            << "# Running Coupling: (RESUM_RC): " << config::RESUM_RC << " (0 fc, 1 balitsky, 2 parent, 3 smallest, 4 guillaume)" << endl
            << "# Running Coupling: (RC_DIS):   " << nlodis_config::RC_DIS << " (0 fc, 1 parent, 2 smallest, 3 guillaume)" << endl
            << "# Use NLOimpact: " << useNLO << endl
            << "# Use SUBscheme: " << useSUB << endl
            << "# Use Sigma3: " << useSigma3 << endl
            << "# Use improved Z2 bound: " << useImprovedZ2Bound << endl
            << "# Use Z2 loop term: " << useBoundLoop << endl
            << "# Cuba MC: " << cubaMethod
                << ", Cuba eps = " << nlodis_config::CUBA_EPSREL
                << ", Cuba maxeval = " << (float)nlodis_config::CUBA_MAXEVAL
                << endl
            << "# config::INTACCURACY = " << config::INTACCURACY
                << ", config::RPOINTS = " << config::RPOINTS
                << ", config::DE_SOLVER_STEP = " << config::DE_SOLVER_STEP
                << ", config::{MINR, MAXR}, nlodis_config::{MINR, MAXR} = " << config::MINR << " " << config::MAXR << " " << nlodis_config::MINR << " " << nlodis_config::MAXR
                << endl;


    return 0;
}
