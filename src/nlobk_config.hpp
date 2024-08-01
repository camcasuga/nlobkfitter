/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

// Configs

#ifndef _NLOBK_CONFIG_HPP
#define _NLOBK_CONFIG_HPP

#include <string>
#include <sstream>



#define LINEINFO __FILE__ << ":" << __LINE__
inline double SQR(double x) { return x*x; }

namespace config
{
    extern bool VERBOSE;
    
    extern double NC;
    extern double NF;

    

    extern double LAMBDAQCD;

    extern int RINTPOINTS;
    extern int THETAINTPOINTS;
    extern double INTACCURACY;
    extern double MCINTACCURACY;
    extern double MAXR;
    extern double MINR;
    extern unsigned int RPOINTS;

    extern size_t MCINTPOINTS;

    extern double DE_SOLVER_STEP;

    // Alpha_s in LO part
    enum RunningCouplingLO
    {
        FIXED_LO,
        PARENT_LO,      // Parent dipole where all beta terms are included
        PARENT_BETA_LO, // Parent dipole, only renormalization scale term is included in as, the second beta term is explicit in the expression
        SMALLEST_LO,
        BALITSKY_LO,
		FRAC_LO, // fastest apparent convergence in 1507.03651
		GUILLAUME_LO // 1708.06557
    };
    enum RunningCouplingNLO
    {
        FIXED_NLO,
        PARENT_NLO,
        SMALLEST_NLO
    };
    extern double FIXED_AS;


    extern RunningCouplingLO RC_LO;
    extern RunningCouplingNLO RC_NLO;

    extern bool DOUBLELOG_LO_KERNEL; // include double log term from the LO kernel
    extern bool ONLY_DOUBLELOG;     // only include double log term

    enum INTEGRATION_METHOD
    {
        VEGAS,
        MISER,
        MULTIPLE            // No monte carlo
    };
    extern INTEGRATION_METHOD INTMETHOD_NLO;

    extern bool ONLY_NLO;   // do not keep as^1 terms

    extern bool FORCE_POSITIVE_N;   // Force N(r)>=0


    extern bool DNDY;   // Print only dn/dy and exit

    extern bool RESUM_DLOG; // Resum double log
    extern bool RESUM_SINGLE_LOG;
    
    enum SINGLELOG_RESUM_RC
    {
		RESUM_RC_FIXED,
		RESUM_RC_BALITSKY,
		RESUM_RC_PARENT,
        RESUM_RC_SMALLEST,
		RESUM_RC_GUILLAUME,
	};
	
	extern SINGLELOG_RESUM_RC RESUM_RC;

    extern bool NO_K2;      // Do not include K_2

    
    extern double KSUB;	// Constant factor in the sigle log resummation log
    
    enum KINEMATICAL_CONSTRAINTS
    {
        KC_BEUF_K_PLUS,
        KC_EDMOND_K_MINUS,
        KC_NONE
    };
    
    extern bool TARGET_KINEMATICAL_CONSTRAINT; // gluon k- < target P-
    
    extern KINEMATICAL_CONSTRAINTS KINEMATICAL_CONSTRAINT; // Solve nonlocal kinematically constrained BK (LO part)
    
    extern bool EULER_METHOD;    // Use Euler method instead of Runge Kutta, must be true if KINEMATICA_CONSTRAINT is used

    const bool LOG_INTERPOLATOR = true; // Flag to determine if interpolate the dipole in log(r), log(N) 
    
}
std::string NLOBK_CONFIG_STRING();

#endif
