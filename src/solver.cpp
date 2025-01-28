/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014
 */

#include "solver.hpp"
#include "dipole.hpp"

#include "nlobk_config.hpp"

#include <cmath>
#include <ctime>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_errno.h>

// Integration constants
const double eps = 1e-30;
using namespace config;
//using std::isinf;
//using std::isnan;
using std::abs;



BKSolver::BKSolver(Dipole* d)
{
    dipole=d;
    tmp_output = "";
    icx0_nlo_impfac=1;
    icTypicalPartonVirtualityQ0sqr=1;
}

struct DEHelper{
    BKSolver* solver;
};

BKSolver::BKSolver()
{
    alphas_scaling=1.0;
    icx0_nlo_impfac=1;
    icTypicalPartonVirtualityQ0sqr=1;
}


int Evolve(double y, const double amplitude[], double dydt[], void *params);
int BKSolver::Solve(double maxy)
{
    /*
     * In order to be able to use Runge Kutta we create one large array
     * which is basically just vector that we evolve
     * Array size is Dipole->RPoints()
     */
    
    cout <<"#### Solving BK equation up to y=" << maxy << endl;
    //cout << "# Nc=" << NC << ", Nf=" << NF << " alphas(r=1) = " << Alphas(1) << endl;
    
    // First check that configs make sense
    if (config::RESUM_DLOG and config::KINEMATICAL_CONSTRAINT != config::KC_NONE)
    {
        cerr << "STOP! User wants to use both resummation of dlogs and kinematical constraint!" << endl;
        exit(1);
    }
    
    
    size_t vecsize = dipole->RPoints();
    double *ampvec = new double [vecsize];
    dipole->InitializeInterpolation(0); // Initialize interpolation at y=yvals[0]=0
    for (unsigned int rind=0; rind<vecsize; rind++)
    {
        ampvec[rind] = dipole->N( dipole->RVal(rind));
    }
    
    double y=0; double step = DE_SOLVER_STEP;  // We have always solved up to y
    int yind=0;
    
    // Intialize GSL
    DEHelper help; help.solver=this;
    gsl_odeiv_system sys = {Evolve, NULL, vecsize, &help};
    
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk2; // rkf45 is more accurate
    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, vecsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (0.00001, 1e-5);    //abserr relerr  
    // Note: absolute uncertainty tolerance here is quite large for the small-r region, but it is in practice needed
    // if the initial condition has a large anomalous dimension: in that case there seems to be quite a fast evolution
    // especially in the region where the dipole is approx. 0
    // Use should check that final results are not sensitive to this parameter!
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (vecsize);
    double h = step;  // Initial ODE solver step size
    
    do
    {
        if (!EULER_METHOD)
        {
            double  nexty = y+step;
            while (y<nexty)
            {
                int status = gsl_odeiv_evolve_apply(e, c, s, &sys,
                                                    &y, nexty, &h, ampvec);
                if (status != GSL_SUCCESS) {
                    cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                    << ": " << gsl_strerror(status) << " (" << status << ")"
                    << " y=" << y << ", h=" << h << endl;
                }
                //if (std::abs(y - (int)(y+0.5))<0.01)
                if(VERBOSE)
                {
                    cout << "\r                                                   " << std::flush;
                    cout << "\r" << "# Evolved up to " << y << "/" << maxy << ", h=" << h << std::flush;
                }
            }
            
            // Check ampvec
            for (int i=0; i<vecsize; i++)
            {
                if (std::isinf(ampvec[i]) or std::isnan(ampvec[i]))
                {
                    cerr << "Ampvec[i=" << i<< "]=" << ampvec[i] << " " << LINEINFO << endl;
                    exit(1);
                }
            }
        }
        else
        {
            // Own implementation of the euler method
            double *dydt = new double[vecsize];
            
            Evolve(y, ampvec, dydt, &help);
            for (int i=0; i<vecsize; i++)
            {
                ampvec[i] = ampvec[i] + step * dydt[i];
            }
            y = y + step;
            delete[] dydt;
            // if(VERBOSE) cout << "# Evolved at y=" << y << endl;
            if(VERBOSE){
                cout << "\r                                               " << std::flush;
                cout << "\r" << "# Evolved at y=" << y << "/" << maxy << std::flush;
            }
        }
        
        yind = dipole->AddRapidity(y, ampvec);
        
        if (tmp_output != "")
            dipole->Save(tmp_output);
        
        // Change Dipole interpolator to the new rapidity
        dipole->InitializeInterpolation(yind);
        
        
    } while (y < maxy);
    
    if(VERBOSE) cout << endl;
    
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    delete[] ampvec;
    return 0;
}

int Evolve(double y, const double amplitude[], double dydt[], void *params)
{
    DEHelper* par = reinterpret_cast<DEHelper*>(params);
    Dipole* dipole = par->solver->GetDipole();
    //cout << "#Evolving, rapidity " << y << endl;
    
    
    // Create interpolators for N(r) and S(r)=1-N(r)
    std::vector<double> rvals,nvals,yvals_s;
    double maxr_interp=-1; // at r>maxr N==1 and S==0
    for (int i=0; i<dipole->RPoints(); i++)
    {
        rvals.push_back(dipole->RVal(i));
        
        double n = amplitude[i];
        if (n>1.0) n=1.0;
        if (n<0 and config::FORCE_POSITIVE_N) n=0;
        nvals.push_back(n);
        
        double s = 1.0-amplitude[i];
        if (s<0) s=0;
        if (s>1 and config::FORCE_POSITIVE_N) s=1.0;
        yvals_s.push_back(s);
    }
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned int i=0; i< dipole->RPoints(); i+=1)
    {
        // It seems to be much more efficeint to initialize interpolators locally
        // for each thread
        Interpolator interp(rvals,nvals,LOG_INTERPOLATOR);
        interp.Initialize();
        interp.SetFreeze(true);
        interp.SetUnderflow(0);
        interp.SetOverflow(1.0);
        //interp.SetMaxX(maxr_interp);
        
        
        
        //if (dipole->RVal(i) < 0.001)
        //    continue;
        // Freeze evolution deep in the saturation region where we know that nothing hapens
        if (amplitude[i] > 0.99999)
         {
             dydt[i]=0;
            continue;
         }
         
        double lo = par->solver->RapidityDerivative_lo(dipole->RVal(i), &interp, y);
        
        double nlo=0;
        if (!NO_K2)
        {
            Interpolator interp_s(rvals,yvals_s, LOG_INTERPOLATOR);
            interp_s.Initialize();
            interp_s.SetFreeze(true);
            interp_s.SetUnderflow(1.0);
            interp_s.SetOverflow(0.0);
            //interp_s.SetMaxX(maxr_interp);
            nlo = par->solver->RapidityDerivative_nlo(dipole->RVal(i), &interp, &interp_s);
        }
        
        if (config::DNDY)
        {
#pragma omp critical
            cout << dipole->RVal(i) << " " << lo << " " << nlo << " " << amplitude[i] << endl;
        }
        dydt[i]= lo + nlo;
        if (std::isnan(dydt[i]) or std::isinf(dydt[i]))
        {
            cerr << "Result " << dydt[i] << " at r " << dipole->RVal(i) << endl;
            dydt[i]=0;
        }
        
    }
    if (config::DNDY)
        exit(1);
    return GSL_SUCCESS;
}



/*
 * *****************
 * LO PART
 * *****************
 */

/*
 * Rapidity derivatives
 * compute \partial_y N(r)
 * First leading order contribution (only one 2d integral), then nlo correction
 *
 * Note: here we assume that dipole amplitude is initialized at correct
 * rapidity so that it can just be evaluated
 */

struct Inthelper_nlobk
{
    BKSolver* solver;
    double r;   // = x - y = parent dipole
    double z;   // = y - z = daughter dipole
    double theta_z; // direction of v
    double z2;   // = y - z' = daughter dipole 2
    double theta_z2; // direction of w
    Interpolator* dipole_interp;
    Interpolator* dipole_interp_s;  // interpolates S=1-N
    double rapidity; // Current rapidity
};

double Inthelperf_lo_z(double v, void* p);
double Inthelperf_lo_theta(double theta, void* p);

// Last argument is optional, and used only with kinematical constraint
double BKSolver::RapidityDerivative_lo(double r, Interpolator* dipole_interp, double rapidity)
{
    gsl_function fun;
    Inthelper_nlobk helper;
    helper.solver=this;
    helper.r=r;
    helper.dipole_interp = dipole_interp;
    helper.rapidity = rapidity;
    
    fun.params = &helper;
    fun.function = Inthelperf_lo_z;
    
    gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(RINTPOINTS);
    
    double minlnr = std::log( 0.5*dipole->MinR() );
    double maxlnr = std::log( 2.0*dipole->MaxR() );
    
    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
                               maxlnr, 0, INTACCURACY, RINTPOINTS,
                               GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status==GSL_ESING)
    {
#pragma omp critical
        cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
        << std::abs(abserr/result) << endl;
    }
    
    
    return result;
}

double Inthelperf_lo_z(double z, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->z=std::exp(z);
    //cout << " r " << helper->r << " v " << helper->v << endl;
    gsl_function fun;
    fun.function=Inthelperf_lo_theta;
    fun.params = helper;
    
    gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(THETAINTPOINTS);
    
    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
                               M_PI, 0, INTACCURACY, THETAINTPOINTS,
                               GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status == GSL_ESING and std::abs(result)>1e-7)
    {
#pragma omp critical
        cerr << "thetaInt failed, z=" << z <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        << std::abs(abserr/result) << endl;
    }
    
    result *= std::exp(2.0*z);  // Jacobian v^2 dv
    result *= 2.0;  // As integration limits are just [0,pi]
    return result;
}

// Rapidity shift Delta_i as in
// http://www.int.washington.edu/talks/WorkShops/int_18_3/People/Triantafyllopoulos_D/Triantafyllopoulos.pdf
// slide "Numerical solution"
// Or new: https://arxiv.org/pdf/1902.06637.pdf eq (5.7)
double RapidityShift(double r, double X)
{
    /*    double A = 1.0;
     double B = 1.0;
     if (r==0)
     return 0;
     
     double kappa = X/r;
     return (1.0 - kappa*kappa) / (A * std::pow(kappa,4) + B*std::pow(kappa,2) + 1.0) * std::log(1.0/(kappa*kappa));
     */
    
    if (X < 1e-10) // Divergence
        return 0; // Would give -inf, force to 0
    
    double result = std::max(0.0, std::log(r*r/(X*X)));
    
    return std::max(0.0, std::log(r*r/(X*X)));
    
}

double StepFunction(double x)
{
    if (x>=0)
        return 1.0;
    else
        return 0.0;
}


double Inthelperf_lo_theta(double theta, void* p)
{
    //cout << " theta " << theta << endl;
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    double r = helper->r;   // x - y
    double z = helper->z;
    
    // X = x - z = r - z
    double Xsqr = r*r + z*z - 2.0*r*z*std::cos(theta);
    
    if (Xsqr < SQR(config::MINR) or z < config::MINR or r < config::MINR)
        return 0;
    
    // we should have Xsqr \approx 0
    double X = std::sqrt( Xsqr );
    // Y = y - z = z
    double Y = z;
    
    
    // Target kinematical constraint
    // Note: if initial condition does not refer to x=1 but some smaller x, this needs to be
    // corrected here when considering the actual gluon p^-
    // This shift we call Y0
    if (config::TARGET_KINEMATICAL_CONSTRAINT)
    {
        double Delta = std::log(1.0/(helper->solver->GetICTypicalPartonVirtualityQ0sqr() * std::min(X*X, Y*Y) + eps));
        
        
        if ( helper->solver->GetX0() * std::exp(-(helper->rapidity - Delta)) > helper->solver->GetICX0_nlo_impfac()) return 0;
    }
    
    
    
    if (std::isnan(X) or std::isnan(Y))
    {
        cerr << "NaN! X=" << X <<", Y=" << Y << " " << LINEINFO << endl;
        exit(1);
    }
    
    if (KINEMATICAL_CONSTRAINT == config::KC_NONE)
    {
        double N_X = helper->dipole_interp->Evaluate(X);    // N(X)
        double N_Y = helper->dipole_interp->Evaluate(Y);
        double N_r = helper->dipole_interp->Evaluate(r);
        
        return helper->solver->Kernel_lo(r, z, theta) * ( N_X + N_Y - N_r - N_X*N_Y );
        
        
    }
    else
    {
        if (!config::EULER_METHOD)
        {
            cerr << "Using KinematicalConstraint but not EulerMethod? " << LINEINFO << endl;
            exit(1);
        }
        
        // Implement kinematical constraint from 1708.06557 Eq. 165
        if (KINEMATICAL_CONSTRAINT == config::KC_BEUF_K_PLUS)
        {
            double delta012 = std::max(0.0, std::log( std::min(X*X, Y*Y) / (r*r) ) ); // (166)
            double shifted_rapidity = helper->rapidity - delta012;
            
            // Step function should respect kinematical boundary x < 1
            if (helper->solver->GetX0() * std::exp(-shifted_rapidity) > helper->solver->GetICX0_nlo_impfac())// Step function in (165)
                return 0;
      
            // Now we are in rapidity "below" initial condition, but still in the kinematically allowed region
            // As we freeze N(r, y < Y_0) = N(r, ic), we can just forget the rapidity shift here
            if (shifted_rapidity < 0)
                shifted_rapidity = 0;
            
            double N_r = helper->dipole_interp->Evaluate(r);
            
            // Dipoles at shifter rapidity
            double s02 = 1.0 - helper->solver->GetDipole()->InterpolateN(X, shifted_rapidity);
            double s12 = 1.0 - helper->solver->GetDipole()->InterpolateN(Y, shifted_rapidity);
            double s01 = 1.0 - N_r;
            
            return helper->solver->Kernel_lo(r, z, theta) * ( -s02*s12 + s01);
        }
        else if (KINEMATICAL_CONSTRAINT == config::KC_EDMOND_K_MINUS)
        {
            
            // Triantafyllopoulos et al, new resummation
            // https://arxiv.org/pdf/1902.06637.pdf
            
            double shifted_1 = helper->rapidity - RapidityShift(r,X);
            double shifted_2 = helper->rapidity - RapidityShift(r,Y);
            
            double icx0 = helper->solver->GetX0();
            double x0 = helper->solver->GetICX0_nlo_impfac();
            double eta0 = std::log(x0/icx0);
            
            // Step function
            // not in (9.3)
            // Note that our evolution starts at eta=0, and eta=0 corresponds to helper->solver->GetX0();
            // So the kinematical requirement for the step function is that xbj < 1, meaning that
            // helper->rapidity - rpaidity_shift > -eta0
            // or helper->rapidity + eta0 - rapidity_shift > 0
            // So note that notation is different than in 1902.06637
            double xyzshift = std::max(0.0, std::log(r*r / std::min( X*X+eps, Y*Y+eps ) ) );
            if (helper->rapidity - xyzshift < 0) return 0;
            
            // (9.3)
            // If shifted rapidity < 0, use initial condition
            double shifted_S_X = 0;
            if (helper->rapidity - RapidityShift(r,X) > 0)
                shifted_S_X = 1.0 - helper->solver->GetDipole()->InterpolateN(X, helper->rapidity - RapidityShift(r,X));
            else
                shifted_S_X = 1.0 - helper->solver->GetDipole()->GetInitialCondition()->DipoleAmplitude(X);
            //shifted_S_X =  1.0 - helper->solver->GetDipole()->InterpolateN(X, helper->rapidity - RapidityShift(r,X));
            
            double shifted_S_Y = 0;
            if (helper->rapidity - RapidityShift(r,Y) > 0)
                shifted_S_Y = 1.0 - helper->solver->GetDipole()->InterpolateN(Y, helper->rapidity - RapidityShift(r,Y));
            else
                shifted_S_Y = 1.0 - helper->solver->GetDipole()->GetInitialCondition()->DipoleAmplitude(Y);
            //        shifted_S_Y = 1.0 - helper->solver->GetDipole()->InterpolateN(Y, helper->rapidity - RapidityShift(r,Y));
            
            
            double S_r = 0;
            if (helper->rapidity > 0)
                S_r = 1.0 - helper->dipole_interp->Evaluate(r);
            else
                S_r = 1.0 - helper->solver->GetDipole()->GetInitialCondition()->DipoleAmplitude(r);
            
            // Check possible numerical errors
            if (shifted_S_X < 0) shifted_S_X  = 0;
            if (shifted_S_Y < 0) shifted_S_Y = 0;
            if (S_r < 0) S_r= 0;
            
            // - as we evolve N, and this is written otherwise of S
            double res =  -helper->solver->Kernel_lo(r, z, theta)  * ( shifted_S_X * shifted_S_Y - S_r );
            if (std::isnan(res))
            {
                cout << "NaN! rapidity " << helper->rapidity << " Xshift " << RapidityShift(r,X) << " Yshift " << RapidityShift(r,Y) << " X=" << X << ", Y=" << Y <<", r=" << r <<", S_X " << shifted_S_X << " S_Y " << shifted_S_Y << endl;
                exit(1);
            }
            return res;
        } // end if KC_EDMOND_K_MINUS
        else
        {
            cerr << "Unknown kinematical constraint specified! " << LINEINFO << endl;
            exit(1);
        }
    }
}


/*
 * LO BK Kernel
 * Evaluates the kernel at given dipole size
 *
 * Parameters:
 *  Parent dipole size r
 *  Daughter dipole size z
 *  Daughter dipole angle theta [0,2\pi]
 */

double BKSolver::Kernel_lo(double r, double z, double theta)
{
    
    // Y = y-z = z
    double Y=z;
    // X = x-z = r - z
    double X = std::sqrt( r*r + z*z - 2.0*r*z*std::cos(theta) );
    
    
    double result=0;
    
    // Handle divergences
    //    if (X<eps or Y<eps)
    //        return 0;
    
    
    // ************************************************** QCD
    
    double min = std::min(X, Y);
    if (r < min)
        min=r;
    
    double alphas_scale = 0;
    
    double eps = 1e-50;
    
    // Note: we have previously checked that X and Y are not zero

    // Fixed as or Balitsky
    // Note: in the limit alphas(r)=const Balitsky -> Fixed coupling as
    if (RC_LO == FIXED_LO )
    {
        result =NC/(2.0*SQR(M_PI))*config::FIXED_AS * SQR(r / (X*Y) );
        alphas_scale=r;
    }
    else if (RC_LO == BALITSKY_LO  )
    {
        double alphas_y = Alphas(Y);
        double alphas_x = Alphas(X);
        result =
        NC/(2.0*SQR(M_PI))*Alphas(r)
        * (
           SQR(r) / ( SQR(X) * SQR(Y)  )
           + 1.0/SQR(Y)*(alphas_y/alphas_x - 1.0)
           + 1.0/SQR(X)*(alphas_x/alphas_y - 1.0)
           );
        alphas_scale = r;
    }
    else if (RC_LO == SMALLEST_LO)
    {
        result = NC*Alphas(min) / (2.0*SQR(M_PI))*SQR(r/(X*Y));
        alphas_scale = min;
    }
    else if (RC_LO == PARENT_LO)
    {
        result = NC*Alphas(r) / (2.0*SQR(M_PI)) * SQR(r/(X*Y));
        alphas_scale = r;
    }
    else if (RC_LO == FRAC_LO)
    {
        // 1507.03651, fastest apparent convergence
        double asbar_r = Alphas(r)*NC/M_PI;
        double asbar_x = Alphas(X)*NC/M_PI;
        double asbar_y = Alphas(Y)*NC/M_PI;
        result = 1.0/(2.0*M_PI) * std::pow(
                                           1.0/asbar_r + (SQR(X)-SQR(Y))/SQR(r) * (asbar_x - asbar_y)/(asbar_x * asbar_y)
                                           , -1.0);
        result = result * SQR(r / (X*Y));
        alphas_scale=r;    // this only affects K1_fin
    }
    else if (RC_LO == GUILLAUME_LO)
    {
        // 1708.06557 Eq. 169
        double r_eff_sqr = r*r * std::pow( Y*Y / (X*X), (X*X-Y*Y)/(r*r) );
        result = NC*Alphas(std::sqrt(r_eff_sqr)) / (2.0*SQR(M_PI)) * SQR(r/(X*Y ));
        alphas_scale = std::sqrt(r_eff_sqr);
        
    }
    else
    {
        cerr << "Unknown LO kernel RC! " << LINEINFO << endl;
        return -1;
    }
    
    if (std::isnan(result) or std::isinf(result))
    {
        //cerr << "Result " << result << " at r=" << r << ", z=" << z << endl;
        result=0;
    }
    
    
    ////// Resummations
    double resummation_alphas = 0;
    if (config::RESUM_RC == RESUM_RC_PARENT)
        resummation_alphas = Alphas(r);
    else if (config::RESUM_RC == config::RESUM_RC_SMALLEST)
        resummation_alphas = Alphas(min);
    else if (config::RESUM_RC == config::RESUM_RC_GUILLAUME)
        resummation_alphas = Alphas(alphas_scale);
    else if (config::RESUM_RC == config::RESUM_RC_BALITSKY)
        cerr << "Check balitsky prescription resummation code! " << LINEINFO << endl;
    else if (config::RESUM_RC == config::RESUM_RC_FIXED)
        resummation_alphas = Alphas(r);
    else
    {
        cerr << "Unknown resummation alphas scale! " << LINEINFO << endl;
        exit(1);
    }
    
    double dlog = 1.0;
    if (config::DOUBLELOG_LO_KERNEL == false or config::RESUM_DLOG == true or config::KINEMATICAL_CONSTRAINT == config::KC_BEUF_K_PLUS)
        dlog=0.0;
    
    
    double resum=1.0;
    if (config::RESUM_DLOG and r > 1.01*config::MINR)
    {
        double x =  4.0*std::log(X/r) * std::log(Y/r) ; // rho^2 in Ref.
        // double x = 2.0*std::log(X/Y); //https://indico.ectstar.eu/event/12/contributions/350/attachments/189/233/2018_ECT_Triantafyllopoulos.pdf
        if (x >=0)
        {
            // argument to the Bessel function is 2sqrt(bar as * x)
            double as_x = std::sqrt( resummation_alphas*NC/M_PI * x );
            resum = gsl_sf_bessel_J1(2.0*as_x) / as_x;
        }
        else // L_xzr L_yzr < 0
        {
            x = std::abs(x);
            double as_x = std::sqrt( resummation_alphas*NC/M_PI * x);
            gsl_sf_result res;
            int status = gsl_sf_bessel_I1_e(2.0*as_x, &res);
            
            if (status != GSL_SUCCESS)
            {
                if (std::isnan(X) or std::isnan(Y))
                    return 0;    // 0/0, probably z=x or z=y
                cerr << "GSL error " << status <<", result " << res.val << ", as_x=" << as_x << ", x=" << x << ", r=" << r<<", X=" << X << ", Y=" << Y << ": " << " z: " << z  << LINEINFO << endl;
                return 0;
            }
            resum = res.val / as_x;
        }
        
        if (std::isnan(resum))
        {
            resum =  1; //1.0;    // 0/0 -> 1 TODO: check
        }
    }
    
    
    
    // Resum single logs
    double singlelog_resum = 1.0;
    double singlelog_resum_expansion = 0;
    double minxy = std::min(X,Y);
    if (std::abs(minxy) < eps) minxy = eps;
    
    if (config::RESUM_SINGLE_LOG)
    {

        double alphabar = resummation_alphas*NC/M_PI;
        const double A1 = 11.0/12.0;
        //singlelog_resum = std::pow( config::KSUB * SQR(r/minxy), sign*A1*alphabar );
        singlelog_resum = std::exp( - alphabar * A1 * std::abs( std::log( config::KSUB * SQR(r/minxy) ) ) );
        
        // remove as^2 part of the single log resummation
        // as it is part of the full NLO coming from K2
        //singlelog_resum_expansion = sign * alphabar * A1 * 2.0*std::log( std::sqrt(config::KSUB)*r/minxy );
        singlelog_resum_expansion = - alphabar * A1 * std::abs( 2.0 * std::log( std::sqrt(config::KSUB) * r/minxy ) ) ;
    }
    
    
    
    if (isnan(result) or isinf(result))
    {
        cerr << "Note: Kernel_lo()=" << result << ", with r=" << r <<", X=" << X << ", Y=" << Y << ", returning 0..." << endl;
        return 0;
    }

    double lo_kernel = Alphas(alphas_scale)*NC/(2.0*M_PI*M_PI) * SQR( r / (X*Y)); // lo kernel with parent/smallest dipole
    double k1fin = lo_kernel * Alphas(alphas_scale) * NC / (4.0*M_PI)
    * (
       67.0/9.0 - SQR(M_PI)/3.0 - 10.0/9.0 * NF/NC
       - dlog*2.0 * 2.0*std::log( X/r ) * 2.0*std::log( Y/r )
       );

    if (KINEMATICAL_CONSTRAINT == config::KC_BEUF_K_PLUS and config::NO_K2 == false) // KCBK + NLO corrections to BK
    {
        if (RESUM_DLOG == true or RESUM_SINGLE_LOG == true)
        {
            cerr << "Kinematical constraint should not be use with double log resummation, single log resummation might be ok... " << LINEINFO << endl;
            exit(1);
        }
        return result + k1fin; // No subtraction term as we don't include the single log resummation
    }
    
    
    if (RESUM_DLOG == false and RESUM_SINGLE_LOG==false) 
    {
        return result;
    }
    if (NO_K2 and (RESUM_DLOG or RESUM_SINGLE_LOG))
    {
        // Resummed K_1, no subtraction or other as^2 terms in K_1
        return resum*singlelog_resum*result;
    }
    
    
    double subtract = 0;
    if (config::RESUM_RC != RESUM_RC_BALITSKY)
        subtract = lo_kernel * singlelog_resum_expansion;
    else    // Balitsky
        subtract = result * singlelog_resum_expansion;
    
    
    
    result = resum*singlelog_resum*result
        - subtract   // remove as^2 part of single log resummation
        + k1fin;
    
    
    return result;
    
    
    
    
    
    
    if (std::isnan(result) or std::isinf(result))
    {
        cerr << "infnan " << LINEINFO << ", r=" << r << ", X=" << X << ", Y=" << Y << endl;
        exit(1);
    }
    
    
    return result;
}




/*
 * *****************
 * NLO PART
 * *****************
 */

double Inthelperf_nlo_z(double v, void* p);
double Inthelperf_nlo_theta_z(double theta, void* p);
double Inthelperf_nlo_z2(double v, void* p);
double Inthelperf_nlo_theta_z2(double theta, void* p);
double Inthelperf_nlo(double r, double z, double theta_z, double z2, double theta_z2, BKSolver* solver, Interpolator* dipole_interp, Interpolator* dipole_interp_s);
double Inthelperf_nlo_mc(double* vec, size_t dim, void* p);

double BKSolver::RapidityDerivative_nlo(double r, Interpolator* dipole_interp, Interpolator* dipole_interp_s)
{
    
    Inthelper_nlobk helper;
    helper.solver=this;
    helper.r=r;
    helper.dipole_interp = dipole_interp;
    helper.dipole_interp_s = dipole_interp_s;
    
    
    double minlnr = std::log( 0.5*dipole->MinR() );
    double maxlnr = std::log( 2.0*dipole->MaxR() );
    
    int status; double  result, abserr;
    
    if (INTMETHOD_NLO == MULTIPLE)
    {
        gsl_function fun;
        fun.params = &helper;
        fun.function = Inthelperf_nlo_z;
        
        gsl_integration_workspace *workspace
        = gsl_integration_workspace_alloc(RINTPOINTS);
        status=gsl_integration_qag(&fun, minlnr,
                                   maxlnr, 0, INTACCURACY, RINTPOINTS,
                                   GSL_INTEG_GAUSS15, workspace, &result, &abserr);
        gsl_integration_workspace_free(workspace);
        
        if (status)
        {
            //#pragma omp critical
            //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
            //<< std::abs(abserr/result) << endl;
        }
    }
    else
    {
        size_t dim=4;
        gsl_monte_function fun;
        fun.params=&helper;
        fun.f = Inthelperf_nlo_mc;
        fun.dim=dim;
        double min[4] = {minlnr, minlnr, 0, 0 };
        double max[4] = {maxlnr, maxlnr, 2.0*M_PI, 2.0*M_PI };
        const gsl_rng_type *T;
        gsl_rng *rnd;
        
        size_t calls = MCINTPOINTS;
        
        
        gsl_rng_env_setup ();
        
        T = gsl_rng_default;
        rnd = gsl_rng_alloc (T);
        
        time_t timer;
        time(&timer);
        int seconds=difftime(timer, 0);
        gsl_rng_set(rnd, seconds);
        
        
        
        const int maxiter_vegas=3;
        
        if (INTMETHOD_NLO == VEGAS)
        {
            gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim);
            gsl_monte_vegas_integrate (&fun, min, max, dim, calls/5, rnd, s,
                                       &result, &abserr);
            //cout <<"#Warmup result " << result << " error " << abserr << endl;
            double prevres = result;
            int iters=0;
            do
            {
                gsl_monte_vegas_integrate (&fun, min, max, dim, calls, rnd, s,
                                           &result, &abserr);
#pragma omp critical
                cout << "#Result(r=" << r <<") " << result << " err " << abserr << " relchange " << (result-prevres)/prevres << " chi^2 " << gsl_monte_vegas_chisq (s) << endl;
                prevres=result;
                iters++;
            }
            while ((std::abs( abserr/result) > 0.3 or std::abs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5 ) and iters<maxiter_vegas );
            //while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5);
            //#pragma omp critical
            if (iters>=maxiter_vegas)
            {
                cout <<"# Integration failed at r=" << r <<", result->0, bestresult "<< result << " relerr " << abserr/result << " chi^2 "  << gsl_monte_vegas_chisq (s) << endl;
                result=0;
            }
            //else
            //    cout << "Integration finished, r=" << r<< ", result " << result << " relerr " << abserr/result << " chi^2 "  << gsl_monte_vegas_chisq (s) << " (intpoints " << calls << ")" << endl;
            gsl_monte_vegas_free(s);
            
        }
        else if (INTMETHOD_NLO == MISER)
        {
            
            // plain or miser
            //gsl_monte_plain_state *s = gsl_monte_plain_alloc (4);
            gsl_monte_miser_state *s = gsl_monte_miser_alloc (4);
            int iter=0;
            
            do
            {
                iter++;
                if (iter>=2)
                {
                    cerr << "Mcintegral didn't converge in 2 iterations (r=" << r << "), result->0 " << LINEINFO << endl;
                    return 0;
                }
                //gsl_monte_plain_integrate
                gsl_monte_miser_integrate
                (&fun, min, max, 4, calls, rnd, s,
                 &result, &abserr);
                //if (std::abs(abserr/result)>0.2)
                //cerr << "#r=" << r << " misermc integral failed, result " << result << " relerr " << std::abs(abserr/result) << ", again.... (iter " << iter << ")" << endl;
            } while (std::abs(abserr/result)>MCINTACCURACY);
            //gsl_monte_plain_free (s);
            gsl_monte_miser_free(s);
            //cout <<"#Integration finished at r=" << r <<", result " << result << " relerr " << abserr/result << " intpoints " << calls << endl;
            
            gsl_rng_free(rnd);
        }
        
    }
    
    
    
    // used for arxiv version evolution for s
    //result *= -SQR(FIXED_AS*NC) / ( 16.0*M_PI*M_PI*M_PI*M_PI);   // alpha_s^2 Nc^2/(16pi^4), alphas*nc/pi=ALPHABAR_s
    // minus sign as we compute here evolution for S=1-N following ref 0710.4330
    
    return result;
}

double Inthelperf_nlo_z(double z, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->z=std::exp(z);
    gsl_function fun;
    fun.function=Inthelperf_nlo_theta_z;
    fun.params = helper;
    
    gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(THETAINTPOINTS);
    
    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
                               2.0*M_PI, 0, INTACCURACY, THETAINTPOINTS,
                               GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }
    
    result *= std::exp(2.0*z);  // Jacobian v^2 dv
    return result;
}

double Inthelperf_nlo_theta_z(double theta, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->theta_z=theta;
    
    gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(RINTPOINTS);
    
    gsl_function fun;
    fun.function=Inthelperf_nlo_z2;
    fun.params=helper;
    
    double minlnr = std::log( 0.5*helper->solver->GetDipole()->MinR() );
    double maxlnr = std::log( 2.0*helper->solver->GetDipole()->MaxR() );
    
    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
                               maxlnr, 0, INTACCURACY, RINTPOINTS,
                               GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }
    
    return result;
    
    
}

double Inthelperf_nlo_z2(double z2, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    helper->z2=std::exp(z2);
    gsl_function fun;
    fun.function=Inthelperf_nlo_theta_z2;
    fun.params = helper;
    
    gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(THETAINTPOINTS);
    
    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
                               2.0*M_PI, 0, INTACCURACY, THETAINTPOINTS,
                               GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    
    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }
    
    result *= std::exp(2.0*z2);  // Jacobian v^2 dv
    return result;
}

double Inthelperf_nlo_theta_z2(double theta_z2, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    double result = Inthelperf_nlo( helper->r, helper->z, helper->theta_z, helper->z2, theta_z2, helper->solver, helper->dipole_interp, helper->dipole_interp_s);
    return result;
}

double Inthelperf_nlo(double r, double z, double theta_z, double z2, double theta_z2, BKSolver* solver, Interpolator* dipole_interp, Interpolator* dipole_interp_s)
{
    // we choose coordinates s.t. y=0 and x lies on positive x axis
    // X = x-z = -z + r
    double X = std::sqrt(r*r + z*z - 2.0*r*z*std::cos(theta_z));
    // Y = y-z = z
    double Y = z;
    // X' = x-z' = r - z'
    double X2=std::sqrt(r*r + z2*z2 - 2.0*r*z2*std::cos(theta_z2) );
    // Y' = y-z' = -z'
    double Y2=z2;
    // z - z'
    double z_m_z2 = std::sqrt( z*z + z2*z2 - 2.0*z*z2*std::cos(theta_z - theta_z2) );
    
    
    double result=0;

    double k = solver->Kernel_nlo(r,X,Y,X2,Y2,z_m_z2);
    double kswap = solver->Kernel_nlo(r,X2,Y2,X,Y,z_m_z2);
    
    
    
    
    /*
     // Dipole part using N
     double dipole = dipole_interp->Evaluate(z_m_z2)
     - dipole_interp->Evaluate(X)*dipole_interp->Evaluate(z_m_z2)
     - dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y2)
     - dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y2)
     + dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y)
     + dipole_interp->Evaluate(X)*dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y2)
     + dipole_interp->Evaluate(Y2) - dipole_interp->Evaluate(Y); // This is not part of eq. (136) in PRD
     
     double dipole_swap = dipole_interp->Evaluate(z_m_z2)
     - dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(z_m_z2)
     - dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y)
     - dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y)
     + dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y2)
     + dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(z_m_z2)*dipole_interp->Evaluate(Y)
     + dipole_interp->Evaluate(Y) - dipole_interp->Evaluate(Y2);
     */
    
    // Dipole part using S, minus sign as the evolution is written for n, not s
    double dipole = -( dipole_interp_s->Evaluate(X)*dipole_interp_s->Evaluate(z_m_z2)*dipole_interp_s->Evaluate(Y2)
                      - dipole_interp_s->Evaluate(X)*dipole_interp_s->Evaluate(Y)  );
    double dipole_swap = -( dipole_interp_s->Evaluate(X2)*dipole_interp_s->Evaluate(z_m_z2)*dipole_interp_s->Evaluate(Y)
                           - dipole_interp_s->Evaluate(X2)*dipole_interp_s->Evaluate(Y2)  );
    
    //result = k*dipole;
    result = (k*dipole + kswap*dipole_swap)/2.0;
    
    if (NF>0)
    {
        double kernel_f = solver->Kernel_nlo_fermion(r,X,Y,X2,Y2,z_m_z2);
        double kernel_f_swap = solver->Kernel_nlo_fermion(r,X2,Y2,X,Y,z_m_z2);
        
        double dipole_f = dipole_interp_s->Evaluate(Y) * ( dipole_interp_s->Evaluate(X2) - dipole_interp_s->Evaluate(X) );
        double dipole_f_swap = dipole_interp_s->Evaluate(Y2) * ( dipole_interp_s->Evaluate(X) - dipole_interp_s->Evaluate(X2) );
        
        /*
         // Dipole part using N
         double dipole_f = dipole_interp->Evaluate(X) - dipole_interp->Evaluate(X2)
         - dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y) + dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y);
         double dipole_f_swap = dipole_interp->Evaluate(X2) - dipole_interp->Evaluate(X)
         - dipole_interp->Evaluate(X2)*dipole_interp->Evaluate(Y2) + dipole_interp->Evaluate(X)*dipole_interp->Evaluate(Y2);
         */
        result += -(kernel_f*dipole_f + kernel_f_swap * dipole_f_swap)/2.0;     // Minus sign as the evolution is written for S and we solve N
        
    }
    
    
    /////// Coefficients
    
    // If the alphas scale is set by the smallest dipole, multiply the kernel here by as^2 and
    // other relevant factors
    
    if (RC_NLO == FIXED_NLO)
        result *= SQR(FIXED_AS*NC) / (8.0*std::pow(M_PI,4) );
    else if (RC_NLO == PARENT_NLO)
        result *= SQR( solver->Alphas(r) * NC) / (8.0 * std::pow(M_PI, 4) );
    else if (RC_NLO == SMALLEST_NLO)
    {
        double min_size = r;
        if (X < min_size) min_size = X;
        if (Y < min_size) min_size = Y;
        if (X2 < min_size) min_size=X2;
        if (Y2 < min_size) min_size = Y2;
        if (z_m_z2 < min_size) min_size = z_m_z2;
        
        result *= SQR( solver->Alphas(min_size) * NC) / (8.0 * std::pow(M_PI, 4) );
    }
    else
    {
        cerr << "Unknown NLO kernel alphas! " << LINEINFO << endl;
        return -1;
    }
        
    
    
    
    
    if (std::isnan(result) or std::isinf(result))
    {
        return 0;
    }
    
    return result;
    
}

/*
 * GSL wrapper for monte carlo integration
 * vec is: [ln u, ln v, theta_u, theta_v]
 */
double Inthelperf_nlo_mc(double* vec, size_t dim, void* p)
{
    Inthelper_nlobk* helper = reinterpret_cast<Inthelper_nlobk*>(p);
    
    if (dim != 4)
    {
        cerr << "Insane dimension at " << LINEINFO << ": " << dim << endl;
        exit(1);
    }
    
    double integrand = Inthelperf_nlo(helper->r, std::exp(vec[0]), vec[2], std::exp(vec[1]), vec[3], helper->solver, helper->dipole_interp, helper->dipole_interp_s)
    * std::exp(2.0*vec[0]) * std::exp(2.0*vec[1]);  // Jacobian
    
    //cout << std::exp(vec[0]) << " " << std::exp(vec[1]) << " " << integrand << endl;
    
    return integrand;
    
}

/***************************************************
 * Non-conformal kernels
 * Note that as^2 nc^2 / (8pi^4) is taken out from the kernels
 ****************************************************/


/*
 * NLO evolution kernel for non-conformal N
 */
double BKSolver::Kernel_nlo(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    
    double kernel = -2.0/std::pow(z_m_z2,4);
    
    kernel += (
               ( SQR(X*Y2) + SQR(X2*Y) - 4.0*SQR(r*z_m_z2) ) / ( std::pow(z_m_z2,4) * (SQR(X*Y2) - SQR(X2*Y)) )
               + std::pow(r,4) / ( SQR(X*Y2)*( SQR(X*Y2) - SQR(X2*Y) )  )
               + SQR(r) / ( SQR(X*Y2*z_m_z2) )
               ) * 2.0*std::log( X*Y2/(X2*Y) );
    
    if (std::isnan(kernel) or std::isinf(kernel))
    {
        //cerr << "Kernel " << kernel <<", r=" << r <<", X=" << X << ", Y=" << Y <<", X2=" << X2 <<", Y2=" << Y2 <<", z-z2=" << z_m_z2 << endl;
        return 0;
    }
    
    return kernel;
}

double BKSolver::Kernel_nlo_fermion(double r, double X, double Y, double X2, double Y2, double z_m_z2)
{
    double kernel=0;
    
    kernel = 2.0 / std::pow(z_m_z2,4.0);
    
    kernel -= (SQR(X*Y2) + SQR(X2*Y) - SQR(r*z_m_z2) ) / ( std::pow(z_m_z2, 4.0)*(SQR(X*Y2) - SQR(X2*Y)) )
    * 2.0*std::log(X*Y2/(X2*Y));
    
    kernel *= NF/NC;        // Divided by NC, as in the kernel we have as^2 nc nf/(8pi^4), but
    // this kernel is multiplied yb as^2 nc^2/(8pi^4)
    
    if (std::isinf(kernel) or std::isnan(kernel))
        return 0;
    
    
    return kernel;
}



////////////////////////////////////////////////
// Running coupling
double BKSolver::Alphas(double r)
{
    if (RC_LO == FIXED_LO or RESUM_RC == RESUM_RC_FIXED or RC_NLO == FIXED_NLO)
        return config::FIXED_AS;
    
    
    double maxalphas=1.0;
    if (config::NF > 3)
    {
        /* Varying n_f scheme (heavy quarks are included), compute effective Lambda_QCD
         * (such that alphas(r) is continuous), see 1012.4408 sec. 2.2.
         */
        
        double dipolescale = 4.0*alphas_scaling/ (r*r);
        double heavyqmasses[2] = {1.3,4.5};
        double nf;
        if (dipolescale < SQR(heavyqmasses[0]))
            nf=3;
        else if (dipolescale < SQR(heavyqmasses[1]))
            nf=4;
        else
            nf = 5;
        
        double b0 = 11.0 - 2.0/3.0*nf;
        // Now we compute "effective" Lambda by requiring that we get the experimental value for alpha_s
        // at the Z0 mass, alphas(Z0)=0.1184, m(Z0)=91.1876
        /*double a0=0.1184;
         double mz = 91.1876;
         double b5 = 11.0 - 2.0/3.0*5.0; // at Z mass all 5 flavors are active
         double b4 = 11.0 - 2.0/3.0*4.0;
         double b3 = 11.0 - 2.0/3.0*3.0;
         double lambda5 = mz * std::exp(-2.0*M_PI / (a0 * b5) );
         double lambda4 = std::pow( heavyqmasses[1], 1.0 - b5/b4) * std::pow(lambda5, b5/b4);
         double lambda3 = std::pow( heavyqmasses[0], 1.0 - b4/b3) * std::pow(lambda4, b4/b3);
         */
        double lambda3=0.146159;
        double lambda4=0.122944;
        double lambda5=0.0904389;
        
        double lqcd;
        if (nf==5) lqcd=lambda5;
        else if (nf==4) lqcd=lambda4;
        else if (nf==3) lqcd=lambda3;
        else
            cerr << "WTF, nf=" << nf <<" at " << LINEINFO << endl;
        double scalefactor = 4.0*alphas_scaling;
        if (scalefactor/(r*r*lqcd*lqcd) < 1.0) return maxalphas;
        double alpha = 4.0*M_PI/(  b0 * std::log(scalefactor/ (r*r*lqcd*lqcd) ) );
        if (alpha > maxalphas) return maxalphas; //NOTE HERE b0 definition have factor 3
        return alpha;
        /*
         double mu0=2.5;
         double c=0.2;
         return 4.0*M_PI / ( b0 * std::log( std::pow( std::pow(mu0, 2.0/c) + std::pow( dipolescale/(lqcd*lqcd), 1.0/c), c) ) );
         */
    }
    
    double Csqr=alphas_scaling;
    double scalefactor = 4.0*Csqr;
    double rsqr = r*r;
    double lambdaqcd=config::LAMBDAQCD;
    const double alphas_mu0=2.5;    // mu0/lqcd
    const double alphas_freeze_c=0.2;
    
    double b0 = (11.0*config::NC - 2.0*config::NF)/3.0;
    
    return 4.0*M_PI / ( b0 * std::log(
                                      std::pow( std::pow(alphas_mu0, 2.0/alphas_freeze_c) + std::pow(4.0*Csqr/(rsqr*lambdaqcd*lambdaqcd), 1.0/alphas_freeze_c), alphas_freeze_c)
                                      ) );
    
    /*
     double scalefactor = 4.0 * alphas_scaling;
     double rsqr = r*r;
     double lambdaqcd=config::LAMBDAQCD;
     if (scalefactor/(rsqr*lambdaqcd*lambdaqcd) < 1.0) return maxalphas;
     double alpha = 12.0*M_PI/( (11.0*NC-2.0*NF)*std::log(scalefactor/ (rsqr*lambdaqcd*lambdaqcd) ) );
     if (alpha>maxalphas)
     return maxalphas;
     return alpha;
     */
}


Dipole* BKSolver::GetDipole()
{
    return dipole;
}

void BKSolver::SetTmpOutput(std::string fname)
{
    tmp_output = fname;
}
