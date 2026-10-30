/*
 * LCPT NLO DIS fitter
 * Henri Hänninen <henri.j.hanninen@student.jyu.fi>, 2017-
 */

#include "nlodis_config.hpp"


// Default configs
namespace nlodis_config
{
    int CUBA_MAXEVAL=1e6;
    double CUBA_EPSREL=0.01;
    bool VERBOSE = false;
    bool PRINTDATA = false;

    double MAXR=50; // 50 for old MV runs
    double MINR=1e-6;

    bool USE_MASSES = false;

    PerfScheme PERF_MODE = nlodis_config::DISABLED;
    RunningCouplingDIS RC_DIS;
    SubtractionScheme SUB_SCHEME;
    QuarkMassScheme MASS_SCHEME;
    SubSchemeTermKernel SUB_TERM_KERNEL;
    GluonZ2Minimum Z2MINIMUM;
    TargetRapidityBKRhoPresc TRBK_RHO_PRESC = nlodis_config::TRBK_RHO_DISABLED;
}
