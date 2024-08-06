/*
 * nloBK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2015
 */

#include "dipole.hpp"
#include "ic.hpp"
#include "amplitudelib/datafile.hpp"
#include "nlobk_config.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <tools/tools.hpp>

using std::cerr; using std::cout; using std::endl;
using namespace config;

/*
 * Constructor
 * Parameter: Initial condition
 * Pointer to the initial condition is saved and used in Save() method.
 */

Dipole::Dipole(InitialCondition* ic_)
{
    X0=1;
    ic=ic_;
    // Initialize rvals
    double step = std::pow(MAXR/MINR, 1.0/RPOINTS);
    std::vector<double> initial_amplitude;
    for (int i=0; i < RPOINTS; i++)
    {
        double r = MINR * std::pow(step, i);
        rvals.push_back(r);
        initial_amplitude.push_back( ic->DipoleAmplitude(r) );
    }
    amplitude.push_back(initial_amplitude);
    
    yvals.push_back(0);

    // Initialize initial condition interpolator
    dipole_interp=NULL;
    InitializeInterpolation(0);
}

/*
 * Creates 1D interpolator which gives dipole amplitude as a function
 * of dipole size at given rapidity yvals[yind]
 *
 */

int Dipole::InitializeInterpolation(int yind)
{
    if (yind >= yvals.size())
    {
        cerr << "Rapidity index " << yind << " is too large, maxindex is " << yvals.size()-1 << LINEINFO << endl;
        return -1;
    }

    if (dipole_interp != NULL) // free old interpolator
    {
        delete dipole_interp;
    }

    dipole_interp = new Interpolator(rvals, amplitude[yind],LOG_INTERPOLATOR);
    dipole_interp->SetUnderflow(0);
    dipole_interp->SetOverflow(1.0);
    dipole_interp->SetFreeze(true);
    dipole_interp->Initialize();
    interpolator_yind = yind;

    return 0;


}

/*
 * Evaluates the initialized interpolator for the given dipole size
 * Trhows an execption and returns -1 if interpolator is not
 * initialized
 */
double Dipole::N(double r)
{
    if (dipole_interp==NULL)
    {
        return -1;
    }

    double n = dipole_interp->Evaluate(r);
    if (n<0 and config::FORCE_POSITIVE_N) return 0;
    if (n>=1.0) return 1.0;

    return dipole_interp -> Evaluate(r);
}

double Dipole::S(double r)
{
    double s=1.0-N(r);
    if (s<0) return 0;
    if (s>1.0 and config::FORCE_POSITIVE_N) return 1.0;
    return s;
}

/*
 * Interpolate dipole amplitude in both r and y
 *
 * Currently, we interpolate linearly in y, and use spline in r
 * We could also use Interpolator2D class, but it requires at least
 * 4 points also in rapidity, and thus would require more extra care
 * close to the initial condition
 */
double Dipole::InterpolateN(double r, double y)
{
    if (y<0) y=0;   // Use initial condition for negative rapidities
    if (y > yvals[yvals.size()-1])
    {
        cerr << "Rapidity " << y << " is larger than current maximum " << yvals[yvals.size()-1] << " " << LINEINFO << endl;
        exit(1);
    }
    
    int yind = FindIndex(y, yvals);
    // FindIndex finds largest index such that yvals[yind] <= y
    if (yind == yvals.size()-1 and std::abs(y - yvals[yvals.size()-1]) < 0.0001)
    {
        // We are asking dipole at the current latest rapidity
        if (interpolator_yind == yind)
            return N(r);
        else
        {
            cerr << "InterpolatorN used to interpolate at the current rapidity, really?? " << LINEINFO << endl;
            exit(1);
        }
    }
    
    // Bilinear interpolation
    int rind = FindIndex(r, rvals);
    
    if (y < 1e-8 and rind+1 < rvals.size())
        return amplitude[0][rind] + (r - rvals[rind]) / (rvals[rind+1]-rvals[rind]) * (amplitude[0][rind+1] - amplitude[0][rind]);
   
   
    // DEBUG for fast testing: bilinear interpolation
	   int rind2 = rind+1;
	   int yind2 = yind+1;
	   if (rind2 >= rvals.size()) return 1.0;
	   double r1_ = rvals[rind];
	   double r2_ = rvals[rind2];

	   if (yind2 >= yvals.size()) return amplitude[yind][rind] + (r - r1_) / (r2_ - r1_) * (amplitude[yind][rind2] - amplitude[yind][rind]);

        double y1_ = yvals[yind];
        double y2_ = yvals[yind2];
        double bilin =  (1.0/( (r2_ - r1_)*(y2_ - y1_) ))*( (amplitude[yind][rind])*(r2_ - r)*(y2_ - y) + (amplitude[yind][rind2])*(r - r1_)*(y2_ - y) + (amplitude[yind2][rind])*(r2_ - r)*(y - y1_) + (amplitude[yind2][rind2])*(r - r1_)*(y - y1_) );
    
    return bilin;
   
    
    // Now we know that yind and yind+1 are acceptable indeces, so we can interpolate in rapidity linearly
    // Construct interpolators in r at both rapidities
    /*
    std::vector<double> rvals_lower;
    std::vector<double> rvals_upper;
    std::vector<double> nvals_lower;
    std::vector<double> nvals_upper;
    
    rind = FindIndex(r, rvals);
    int min_index = std::max(rind-3, 0);
    if (min_index + 6 > rvals.size())
    {
        // So larger, that it is good approximation to say that N=1
        return 1.0;
    }
    
    for (int i=0; i<6; i++)
    {
        rvals_lower.push_back(rvals[min_index+i]);
        nvals_lower.push_back(amplitude[yind][min_index+i]);
        rvals_upper.push_back(rvals[min_index+i]);
        nvals_upper.push_back(amplitude[yind+1][min_index+i]);
    }
    Interpolator lower(rvals_lower, nvals_lower);
    Interpolator upper(rvals_upper, nvals_upper);
    double n_lower = lower.Evaluate(r);
    double n_upper = upper.Evaluate(r);
    
    // Linear interpolation in y
    double result =  n_lower + (y - yvals[yind]) / ( yvals[yind+1] - yvals[yind]) * (n_upper - n_lower);
    
    // Possible rounding errors
    if (result > 1.0) return 1.0;
    if (result < 0) return 0;
    
    
    return result;
    */
}

/*
 * Add new rapidity
 * The DE solver in BKSolver class obtains the new dipole amplitude
 * at different dipole sizes and saves it
 * Here we get the table containing the dipole amplitudes and the rapidity value
 *
 * Returns -1 if an error occurs, the corresponding rapidity index if
 * succesfull (that is, YVal[returned_vlaue)=y
 */
int Dipole::AddRapidity(double y, double rgrid[])
{

    if (y<=yvals[yvals.size()-1])
    {
        cerr << "Too small rapidity " << y << " at " << LINEINFO << endl;
        return -1;
    }

    yvals.push_back(y);
    std::vector<double> amp;
    for (unsigned int i=0; i<RPoints(); i++)
    {
        amp.push_back(rgrid[i]);
    }

    amplitude.push_back(amp);

    return yvals.size()-1;

}

/*
 * Save amplitude to the given file
 *
 * Syntax is the same as in rbk and AmplitudeLib
 */
std::string NLOBK_CONFIG_STRING();
int Dipole::Save(std::string filename)
{
    std::ofstream out;
    out.open(filename.c_str());

    // Save info
    out << "# NLO BK equation solver " << VERSION << " build " << " (build " <<  __DATE__ << " " << __TIME__ << ")" << endl;
    out <<"# Initial condition: " << ic->GetString() << endl;

    out << "# " << NLOBK_CONFIG_STRING() << endl;
    
    out << "###" << std::scientific << std::setprecision(15) << MinR() << endl;
    out << "###" << std::scientific << std::setprecision(15) <<
        rvals[1]/rvals[0]  << endl;   // rmultiplier
    out << "###" << RPoints() << endl;
    out << "###" << X0 << endl;   // x0

    for (int yind=0; yind<=YPoints()-1; yind++)
    {
        out << "###" << std::scientific << std::setprecision(15)
            << yvals[yind] << endl;
        for (int rind=0; rind<RPoints(); rind++)
        {
            out << std::scientific << std::setprecision(15)
                << amplitude[yind][rind] << endl;
        }
    }
    out.close();
    return 0;
}


double Dipole::MinR()
{
    return rvals[0];
}
double Dipole::MaxR()
{
    return rvals[rvals.size()-1];
}

unsigned int Dipole::RPoints()
{
    if (rvals.size() != RPOINTS)
    {
        cerr << "STOP! rvals.size() != RPOINTS" << endl;
        exit(1);
    }
    return rvals.size();
}

unsigned int Dipole::YPoints()
{
    return yvals.size();
}

double Dipole::RVal(unsigned int rind)
{
    if (rind < 0 or rind>rvals.size()-1)
    {
        cerr << "Invalid dipole size index " << rind << " at " << LINEINFO << endl;
        return 0;
    }
    return rvals[rind];
}

double Dipole::YVal(unsigned int yind)
{
    return yvals[yind];
}

Dipole::~Dipole()
{
    // if(VERBOSE)     cout << "# Destroying dipole..." << endl;
    if (dipole_interp != NULL)
    {
        delete dipole_interp;
    }
}


Dipole::Dipole(std::string filename)
{
    DataFile data(filename);
    data.GetData(amplitude, yvals);
    int rpoints = data.RPoints();
    double rmultiplier = data.RMultiplier();
    double minr= data.MinR();
    for (int i=0; i<=rpoints; i++)
    {
        rvals.push_back(minr * std::pow(rmultiplier, i));
    }
    if (rvals.size() != amplitude[0].size())
    {
        cerr<< "I got " << rvals.size() << " r values and " <<amplitude[0].size() << " dipoles at y=0!" << LINEINFO << endl;
        exit(1);
    }
    X0=data.X0();
    cout << "# Data read from file " << filename << " x0=" << data.X0() << "  maxy " << yvals[yvals.size()-1] << " rpoints/rapidity " << amplitude[0].size() << endl;
}
