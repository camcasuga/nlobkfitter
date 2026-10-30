
#include <iostream>
#include <amplitudelib/amplitudelib.hpp>
#include <tools/interpolation.hpp>
#include "solver.hpp"
#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"
#include <csignal>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <ctime>
#include <unistd.h>
#include <tools/tools.hpp>
#include <sstream>
#include "nlodissigmar.hpp"

void ErrHandler(const char * reason,
                const char * file,
                int line,
                int gsl_errno);

using namespace std;


int main()
{
    gsl_set_error_handler(&ErrHandler);

    // ***Fit parameters***
    // Initial condition
    double qs0sqr = 0.241*0.241*4;    // Q_s,0^2 at x=0.01 (GeV^2)
    double e_c = 1.0;
    double anomalous_dimension = 1.0;  // probably want to keep constant

    // BK solver
    double alphas_scaling = exp(-2.0*0.5772);     // C^2 in the expression for alpha_s


    // ***Set other configurations, should not be changed during the fit process***
    double maxy = 5.2;      // Solve BK up to this rapidity

    config::RC_LO = config::BALITSKY_LO; // Balitsky running coupling for LO kernel
    config::RESUM_RC = config::RESUM_RC_PARENT; // Parent dipole in the resummation
    config::RESUM_DLOG = true; // Resum doulbe logs
    config::RESUM_SINGLE_LOG = true; // Resum single logs
    config::LO_BK = true;  // Solve LO BK with running coupling, overrides RESUM settings
    config::KSUB = 0.65;  // Optimal value for K_sub
    config::NO_K2 = true;  // Do not include numerically demanding full NLO part
    config::INTACCURACY = 0.02;
    config::MINR = 1e-5;
    config::MAXR = 20;
    config::RPOINTS = 100;

    // Constants
    config::NF= 3;//3;   // Only light quarks
    config::LAMBDAQCD = 0.241;

    // ***Solve resummed BK***
    MV ic;   // Initial condition
    ic.SetQsqr(qs0sqr);
    ic.SetAnomalousDimension(anomalous_dimension);
    ic.SetE(e_c);       // e_c of MVe parametrization

    Dipole dipole(&ic);


    BKSolver solver(&dipole);
    solver.SetAlphasScaling(alphas_scaling);
    solver.Solve(maxy);  // Solve up to maxy


    // ***Study solution***

    // Give solution to the AmplitudeLib object
    AmplitudeLib DipoleAmplitude(solver.GetDipole()->GetData(), solver.GetDipole()->GetYvals(), solver.GetDipole()->GetRvals());

    // Test: print dipole amplitude at rapidities 0,0.5,1

    cout << "# r   N(y=0)  N(y=1)   N(y=5)" << endl;
    for (double r=0.0001; r<10; r*=1.01)
    {
        cout << r<< " " << DipoleAmplitude.N(r, 0.01*exp(-0)) << " " << DipoleAmplitude.N(r, 0.01*exp(-1)) << " " << DipoleAmplitude.N(r,0.01*exp(-5)) << endl;
    }

    SetDipole(&AmplitudeLib::N, DipoleAmplitude);
    SetRunningCoupling(alpha_bar_running_pd);
    std::cout << "SigmarLO(1,1e-3)=" << SigmarLO(1,1e-3,1) << '\n';


    // Test: print saturation scale and its evolution speed
    /*
    vector<double> yvals;
    vector<double> lnqsqrvals;

    for (double y=0; y<maxy-0.1; y+=0.1)
    {
        double xbj =0.01*exp(-y);
        DipoleAmplitude.InitializeInterpolation(xbj);
        double qs = 2.0 / pow(DipoleAmplitude.SaturationScale(xbj, 1.0 - exp(-0.5)),2.0);
        yvals.push_back(y);
        lnqsqrvals.push_back(log(qs));

    }

    Interpolator satscaleinterp(yvals, lnqsqrvals);

    for (double y=0; y<maxy-0.1; y+=0.1)
        cout << y << " " << exp(0.5*satscaleinterp.Evaluate(y)) << " " << satscaleinterp.Derivative(y) << endl;
    */

    return 0;
}
