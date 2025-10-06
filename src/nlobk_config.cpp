/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2015  
 */

#include "nlobk_config.hpp"
#include <sstream>
#include <string>
#include <iostream>

using namespace std;

using namespace config;

// Default configs
namespace config
{
     bool VERBOSE = true;

     double NC=3;
     double NF=3;

     double LAMBDAQCD = 0.241;

     int RINTPOINTS= 85;
     int THETAINTPOINTS = 85;
     double INTACCURACY=0.005;
     double MCINTACCURACY = 0.2;
     double MAXR = 10;          // Quite small, only for testing
     double MINR= 1e-5;//1e-6;
     unsigned int RPOINTS = 170; //170

     size_t MCINTPOINTS = 1e6;

     double DE_SOLVER_STEP = 0.2; // 0.05 paperissa


     double FIXED_AS = 0.2;


     RunningCouplingLO RC_LO = BALITSKY_LO;
     RunningCouplingNLO RC_NLO = PARENT_NLO;

     bool DOUBLELOG_LO_KERNEL = false; // include double log term from the LO kernel
     bool ONLY_DOUBLELOG = false;

     INTEGRATION_METHOD INTMETHOD_NLO = VEGAS;

     bool FORCE_POSITIVE_N = true;

     bool DNDY=false;

     bool ONLY_NLO = false;

     bool RESUM_DLOG = false;
     bool RESUM_SINGLE_LOG = false;

     bool NO_K2 = false;
     
     double KSUB = 1.0;
     
     SINGLELOG_RESUM_RC RESUM_RC = RESUM_RC_PARENT;
    
     KINEMATICAL_CONSTRAINTS KINEMATICAL_CONSTRAINT = KC_NONE;
    
    bool TARGET_KINEMATICAL_CONSTRAINT=false;
    
    bool EULER_METHOD = false;
}


std::string NLOBK_CONFIG_STRING()
{
    std::stringstream ss;
    
    
    
    
    ss<< "# LO Kernel RC: ";
    if (RC_LO == FIXED_LO )
    ss << " fixed as=" << FIXED_AS;
    else if (RC_LO == SMALLEST_LO)
    ss << " smallest dipole";
    else if (RC_LO == BALITSKY_LO)
    ss << " Balitsky";
    else if (RC_LO == PARENT_LO)
    ss << " Parent dipole";
    else if (RC_LO == PARENT_BETA_LO)
    ss << " Parent dipole, explicit beta";
    else if (RC_LO == GUILLAUME_LO)
    ss << " Guillaume RC (arXiv:1708.06557)";
    else if (RC_LO == FRAC_LO)
    ss << " Fastest apparent convergence";
    else
    ss << " NO STRING IMPLEMENTED!";
    
    
    ss<< ". NLO Kernel RC: ";
    if (RC_NLO == FIXED_NLO)
    ss << " fixed as=" << FIXED_AS;
    else if (RC_NLO == SMALLEST_NLO)
    ss << " smallest dipole";
    else if (RC_NLO  == PARENT_NLO)
    ss << " Parent dipole";
    else
    ss << " NO STRING IMPLEMENTED!";
    ss << endl;
    ss <<"# Nc=" << NC << ", Nf=" << NF << endl;
    
    if (config::EULER_METHOD)
        ss <<"# Using Euler method in K1" << endl;
    else
        ss << "# Using RungeKutta method in K1" << endl;
    
    
    //if (FORCE_POSITIVE_N)
    //ss << "# Amplitude is limited to [0,1]" << endl;
    

    //BKSolver sol;
    //ss << "# Alphas(r=1 GeV^-1) = " << sol.Alphas(1) << endl;
    
    if (config::ONLY_NLO) ss << "# keeping only NLO terms" << endl;
    if (config::RESUM_DLOG)
        ss << "# Resumming double log"<< endl;
    
    if (config::RESUM_SINGLE_LOG)
    {
        
        ss << "# Resumming single log, K_sub=" << config::KSUB;
        if (config::RESUM_RC == RESUM_RC_PARENT) ss << " resum rc: parent";
        else if (config::RESUM_RC == RESUM_RC_SMALLEST) ss << " resum rc: smallest";
        else if (config::RESUM_RC == RESUM_RC_BALITSKY) ss << " resum rc: balitsky";
        ss << endl;
    }
    
    ss << "# BK K1 integration relative accuracy: " << INTACCURACY ;
    if (config::NO_K2)
    {
        ss <<  "# Not including K2 and Kf" << endl;
    }
    else
    {
        ss << "# NLO BK as^2 terms: MC integration method: ";
        if (INTMETHOD_NLO == MISER)
            ss <<"Miser, points=" << MCINTPOINTS;
        else if (INTMETHOD_NLO == VEGAS)
            ss <<"Vegas, points=" << MCINTPOINTS;
        else if (INTMETHOD_NLO == MULTIPLE)
            ss << "Multiple integrals (no montecarlo)";
        else
            ss <<"UNKNOWN!";
        ss << endl;
    }
    
    ss <<  "# Kinematical constraint: ";
    if (config::KINEMATICAL_CONSTRAINT == config::KC_NONE)
        ss << "none";
    else if (config::KINEMATICAL_CONSTRAINT == config::KC_BEUF_K_PLUS)
        ss << "Beuf KC (arXiv:1708.06557), evolution in k^+" ;
    else if (config::KINEMATICAL_CONSTRAINT == config::KC_EDMOND_K_MINUS)
        ss << "Edmond KC (arXiv:1902.06637), evolution in k^-";
    else ss << "UNKNOWN!";
    ss << endl;
    
    ss << "# Target kinematical constraint: ";
    if (config::TARGET_KINEMATICAL_CONSTRAINT)
        ss << "enabled" << endl;
    else
        ss << "disabled" << endl;
    
    
    
    return ss.str();
}

