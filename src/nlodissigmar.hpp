
#ifndef _NLODIS_SIGMAR_H
#define _NLODIS_SIGMAR_H

#include <string>
#include <vector>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnUserParameterState.h>

#include <amplitudelib/amplitudelib.hpp>
#include <tools/interpolation.hpp>
#include "solver.hpp"
#include "ic.hpp"
#include "mv.hpp"
#include "ic_datafile.hpp"
#include "dipole.hpp"
#include "solver.hpp"
#include "nlodis_config.hpp"

#include "cuba-4.2.h"
#include "data.hpp"

#define PARALLEL_CHISQR


using namespace ROOT::Minuit2;
using namespace std;

static inline double Sq(double x){return x*x;}

//namespace sigmar_config{extern double maxy;}

class NLODISFitter : public FCNBase
{
public:
    // MINUIT functions
    double operator() (const vector<double>& par) const; // Calculate Chi^2
    double Up() const {return 1.;} // one-standard-deviation errors for fit parameters.

    // Initialize based on MINUIT parameters
    NLODISFitter(MnUserParameters parameters_);

    // Class configurators
    void AddDataset(Data& d);
    void SetNLO(bool s) { computeNLO = s; }
    void SetSUB(bool s) { UseSub = s; }
    void SetSigma3(bool s) { UseSigma3 = s; }
    void UseImprovedZ2Bound(bool b) { useImprovedZ2Bound = b;}
    void UseConsistentlyBoundLoopTerm(bool b) { useBoundLoop = b;}
    void SetCubaMethod(string s) { cubaMethod = s; }


private:

    bool computeNLO, UseSub, UseSigma3, useImprovedZ2Bound, useBoundLoop;
    string cubaMethod;
    MnUserParameters parameters;
    vector<Data*> datasets;

};


class ComputeSigmaR
{
public:
  typedef  double (ComputeSigmaR::*CmptrMemFn)(double x);
  typedef  double (ComputeSigmaR::*CmptrMemFn_void)(void *x);
  typedef  double (ComputeSigmaR::*z2funpointer)(double x, double q);
  typedef  double (ComputeSigmaR::*subterm_z2upper_fp)(double z);
  typedef  double (ComputeSigmaR::*xrapidity_funpointer)(double x, double q);
  typedef  double (ComputeSigmaR::*xrapidity_y_eta_funpointer)(double x, double q, double rsq);
  typedef  double (ComputeSigmaR::*xrapidity_NLO_funpointer)(double Qsq, double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq);
  typedef  double (ComputeSigmaR::*bkkernel_funpointer)(double x01sq, double x02sq, double x21sq);
  typedef  double (ComputeSigmaR::*trbk_rho_funpointer)(double Qsq, double x01sq, double x02sq, double x21sq);
  // old
  typedef  double (ComputeSigmaR::*INLOqg_subterm_fp)(double Q, double x, double z1, double z2,
                                                      double x01sq, double x02sq, double phix0102);


    ComputeSigmaR(AmplitudeLib *ObjectPointer);

    double SigmarLO ( double Q , double xbj, double y) ;
    double SigmarLOmass ( double Q , double xbj, double y, bool charm=false) ;
    double SigmarNLOunsub ( double Q , double xbj, double y) ;
    double SigmarNLOunsub_massive ( double Q , double xbj, double y, double mf) ;
    double SigmarNLOsub ( double Q , double xbj, double y) ;
    double SigmarNLOunsub_UniformZ2Bound ( double Q , double xbj, double y) ;
    double SigmarNLOsub_UniformZ2Bound ( double Q , double xbj, double y) ;
    double SigmarNLOunsub_sigma3 ( double Q , double xbj, double y) ;

    double SigmarNLOunsubRisto ( double Q , double xbj, double y) ;
    double SigmarNLOsubRisto ( double Q , double xbj, double y) ;

    double Structf_LLO ( double Q , double xbj) ;
    double Structf_TLO ( double Q , double xbj) ;
    double Structf_LNLOdip ( double Q , double xbj) ;
    double Structf_TNLOdip ( double Q , double xbj) ;
    double Structf_LNLOdip_z2 ( double Q , double xbj) ;
    double Structf_TNLOdip_z2 ( double Q , double xbj) ;
    double Structf_LNLOqg_unsub ( double Q , double xbj) ;
    double Structf_TNLOqg_unsub ( double Q , double xbj) ;
    double Structf_LNLOqg_sub ( double Q , double xbj) ;
    double Structf_TNLOqg_sub ( double Q , double xbj) ;
    double Structf_LNLOsigma3 ( double Q , double xbj) ;
    double Structf_TNLOsigma3 ( double Q , double xbj) ;

    double Structf_LLO_massive ( double Q , double xbj, double mf) ;
    double Structf_TLO_massive ( double Q , double xbj, double mf) ;
    double Structf_LNLOdip_massive ( double Q , double xbj, double mf) ;
    double Structf_TNLOdip_massive ( double Q , double xbj, double mf) ;
    double Structf_LNLOdip_z2_massive ( double Q , double xbj, double mf) ;
    double Structf_TNLOdip_z2_massive ( double Q , double xbj, double mf) ;
    double Structf_LNLOqg_unsub_massive ( double Q , double xbj, double mf) ;
    double Structf_TNLOqg_unsub_massive ( double Q , double xbj, double mf) ;

    double Structf_LFULLNLOunsub ( double Q , double xbj) ; // no dipole term in these 4 for now.
    double Structf_TFULLNLOunsub ( double Q , double xbj) ;
    double Structf_LFULLNLOsub ( double Q , double xbj) ;
    double Structf_TFULLNLOsub ( double Q , double xbj) ;

    void SetQuarkMassLight(double qMass_light_){ qMass_light = qMass_light_; }
	void SetQuarkMassCharm(double m) { qMass_charm = m; }
    void SetAlphasScalingC2(double c2_){ alpha_scaling_C2_ = c2_; }
    void SetX0(double x0_){ icX0 = x0_; }
    void SetX0_BK(double x0_){ icX0_bk = x0_; }
    void SetQ0Sqr(double q0_){ icQ0sqr = q0_; }
    void SetY0(double y0_){ icY0 = y0_; }

    void SetRunningCoupling(CmptrMemFn p){AlphabarPTR = p;} // function pointer setter
    void SetRunningCoupling_QG(CmptrMemFn_void p){Alphabar_QG_PTR = p;} // function pointer setter

    void SetImprovedZ2Bound(z2funpointer p){z2limit_PTR = p;} // function pointer setter
    void SetEvolutionX_LO(xrapidity_y_eta_funpointer p){Xrpdty_LO_PTR = p;}
    void SetEvolutionX_LO_z2scheme(xrapidity_funpointer p){Xrpdty_LO_projectileY_z2min_PTR = p;}
    void SetEvolutionX_DIP(xrapidity_y_eta_funpointer p){Xrpdty_DIP_PTR = p;}
    void SetEvolutionX_NLO(xrapidity_NLO_funpointer p){Xrpdty_NLO_PTR = p;}

    void SetSigma3BKKernel(bkkernel_funpointer p){K_kernel_PTR = p;} // function pointer setter
    void SetSubTermKernel(nlodis_config::SubSchemeTermKernel sstk){
        if (sstk == nlodis_config::SUBTERM_LOBK_Z2TOZERO){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_lobk_z2tozero;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_lobk_z2tozero;
        }else if (sstk == nlodis_config::SUBTERM_LOBK_EXPLICIT){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_lobk_explicit;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_lobk_explicit;
        }else if (sstk == nlodis_config::SUBTERM_RESUM){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_resumbk;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_resumbk;
        }else if (sstk == nlodis_config::SUBTERM_KCBK_BEUF){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_kcbk_beuf;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_kcbk_beuf;
        }else if (sstk == nlodis_config::SUBTERM_TRBK_EDMOND){
            ILNLOqg_subterm_PTR = &ComputeSigmaR::ILNLOqg_subterm_trbk_edmond;
            ITNLOqg_subterm_PTR = &ComputeSigmaR::ITNLOqg_subterm_trbk_edmond;
        }
    }
    void SetTRBKRhoPrescription(nlodis_config::TargetRapidityBKRhoPresc rhopresc){
        if (rhopresc == nlodis_config::TRBK_RHO_QQ0){
	        rho_PTR = &ComputeSigmaR::rho_rapidity_shift_QQ0;
	        cout << "# TRBK RHO is TRBK_RHO_QQ0" << endl; }
        else if (rhopresc == nlodis_config::TRBK_RHO_RQ0){
	        rho_PTR = &ComputeSigmaR::rho_rapidity_shift_RQ0;
	        cout << "# TRBK RHO is TRBK_RHO_RQ0" << endl; }
        else if (rhopresc == nlodis_config::TRBK_RHO_X_R){
	        rho_PTR = &ComputeSigmaR::rho_rapidity_shift_XR;
	        cout << "# TRBK RHO is TRBK_RHO_X_R" << endl;
	        cout << "############### SHOULD BE WELL CHECKED BEFORE USE ##############" << endl; }
        else if (rhopresc == nlodis_config::TRBK_RHO_MAX_X_Y_R){
	        rho_PTR = &ComputeSigmaR::rho_rapidity_shift_MAX_XYR;
	        cout << "# TRBK RHO is TRBK_RHO_MAX_X_Y_R" << endl;
	        cout << "############### SHOULD BE WELL CHECKED BEFORE USE ##############" << endl; }
        else if (rhopresc == nlodis_config::TRBK_RHO_DISABLED){
            rho_PTR = &ComputeSigmaR::rho_no_shift_projectileY;
	        cout << "# TRBK RHO is TRBK_RHO_DISABLED" << endl;
        }
    }
    void MetaPrescriptionSetter(){
        // NLO: set runningcoupling and C2=Csq for the object.
        ComputeSigmaR::CmptrMemFn alphas_temppointer;
        ComputeSigmaR::CmptrMemFn_void alphas_temppointer_QG;
        if      (nlodis_config::RC_DIS == nlodis_config::DIS_RC_FIXED){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_fixed;
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_fixed;
            cout << "# Using FIXED_LO" << endl;}
        else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_PARENT){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd;
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_pd;
            cout << "# Using parent dipole RC" << endl;}
        else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_GUILLAUME){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd;
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_guillaume;
            cout << "# Using Guillaume RC" << endl;}
        else if (nlodis_config::RC_DIS == nlodis_config::DIS_RC_SMALLEST){
            alphas_temppointer = &ComputeSigmaR::alpha_bar_running_pd; // parent dipole in the LO-like qqbar terms where there are no daughter dipoles.
            alphas_temppointer_QG  = &ComputeSigmaR::alpha_bar_QG_running_smallest;
            cout << "# Using smallest dipole RC" << endl;}
        else {
            cout << "ERROR: Problem with the choice of runnincoupling. Unkonwn nlodis_config::RC_DIS." << endl;
            exit(1);
            }
        this->SetRunningCoupling(alphas_temppointer);
        this->SetRunningCoupling_QG(alphas_temppointer_QG);

        // Set z2 lower limit settings
        ComputeSigmaR::z2funpointer z2bound_funptr;
        if  (nlodis_config::Z2MINIMUM == nlodis_config::Z2IMPROVED) {
            z2bound_funptr = &ComputeSigmaR::z2bound_improved;
            cout << "z2_min scheme is improved z2min(Q^2) ~ 1/Q^2" << endl;
        } else {
            z2bound_funptr = &ComputeSigmaR::z2bound_simple;
            cout << "z2_min scheme is simple z2min = xbj/X_0if" << endl;
        }
        this->SetImprovedZ2Bound(z2bound_funptr);

        // Dipole evalution X, (y = log(x0/X))
        // if using 'sub' scheme with z2imp one must use the extended evolution variable for sigma_LO as well.
        ComputeSigmaR::xrapidity_funpointer x_Y_z2min_ptr;
        if (nlodis_config::SUB_SCHEME == nlodis_config::SUBTRACTED 
            and nlodis_config::Z2MINIMUM == nlodis_config::Z2IMPROVED) {
                x_Y_z2min_ptr = &ComputeSigmaR::Xrpdty_LO_improved;
                cout << "# Using SUB + Z2IMP -- take Q^2 improved z2min scheme into account in LO rapidity!" << endl;
                }
        else {
            x_Y_z2min_ptr = &ComputeSigmaR::Xrpdty_LO_simple;
            cout << "# Using Z2SIM LO rapidity -- simple Q^2 independent z2min scheme in LO rapidity (not using sub && z2imp)" << endl;
            }
        
        ComputeSigmaR::xrapidity_y_eta_funpointer x_lo_y_eta_rap_ptr; // only sub scheme LO term should have kinematical rapidity shift. No shift with unsub IC term!
        ComputeSigmaR::xrapidity_y_eta_funpointer x_dip_y_eta_rap_ptr; //
        ComputeSigmaR::xrapidity_NLO_funpointer x_nlo_fun_ptr;
        if (config::KINEMATICAL_CONSTRAINT == config::KC_EDMOND_K_MINUS) {
            if (nlodis_config::SUB_SCHEME == nlodis_config::SUBTRACTED){
                // sub scheme has evolution in the LO term so it needs the kinematical shift there as well.
                x_lo_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_targetETA;
                x_dip_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_targetETA;
                cout << "# Using shifted target ETA rapidity in SUB leading order term AND NLO dipole term." << endl;
                x_nlo_fun_ptr = &ComputeSigmaR::Xrpdty_NLO_targetETA;
                cout << "# Using target ETA evolution rapidity in NLO terms" << endl;
            }else if (nlodis_config::SUB_SCHEME == nlodis_config::UNSUBTRACTED){
                // unsub scheme has no evolution in the lowest order term, so no shift there, use the Y rapidity technology that doesn't shift.
                // Initial condition will be defined in x_eta = x0_bk.
                x_lo_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_projectileY;
                cout << "# Using unshifted target ETA rapidity in UNSUB lowest order term. (sig^IC no shift)" << endl;
                x_dip_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_targetETA;
                cout << "# Using shifted target ETA rapidity in virtual dipole term " << endl;
                // x_dip_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_projectileY;
                // cout << "# !!!!!! Using UNshifted rapidity in virtual dipole term !!!!!!!!" << endl;
                x_nlo_fun_ptr = &ComputeSigmaR::Xrpdty_NLO_targetETA;
                cout << "# Using target ETA evolution rapidity in NLO terms" << endl;
            }else{
                cout << "# SUB_SCHEME not valid. Exit.";
                exit(0);
            }
        }
        else {
            // Not using Target Rapidity -> no shifts / change of variables to eta
            x_lo_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_projectileY;
            x_dip_y_eta_rap_ptr = &ComputeSigmaR::Xrpdty_LO_projectileY;
            x_nlo_fun_ptr = &ComputeSigmaR::Xrpdty_NLO_projectileY;
            cout << "# Using projectile Y evolution rapidity" << endl;
            // cov_to_eta_x calls rho_shift
            // rho_PTR = &ComputeSigmaR::rho_no_shift_projectileY; // set no_shift here or in SetTRBKRhoPresc
            }
        this->SetEvolutionX_LO(x_lo_y_eta_rap_ptr);
        this->SetEvolutionX_DIP(x_dip_y_eta_rap_ptr);
        this->SetEvolutionX_LO_z2scheme(x_Y_z2min_ptr);
        this->SetEvolutionX_NLO(x_nlo_fun_ptr);

        // sigma_3 BK correction
        this->SetSigma3BKKernel(&ComputeSigmaR::K_resum);

        // sub scheme subtraction term choice
        this->SetSubTermKernel(nlodis_config::SUB_TERM_KERNEL);
        
        // trbk rho prescription
        this->SetTRBKRhoPrescription(nlodis_config::TRBK_RHO_PRESC);
    }

    void SetCubaMethod(string s){cubamethod = s;}

	AmplitudeLib* GetDipole() { return ClassScopeDipolePointer; }

//private:
    //variables
    AmplitudeLib *ClassScopeDipolePointer;
    double qMass_light, alpha_scaling_C2_, icX0, icX0_bk, icY0, icQ0sqr;
	double qMass_charm;
    struct QMasses{
        double m_u, m_d, m_s, m_c, m_b, m_t;
    };
    
    // Method and function pointers
    CmptrMemFn AlphabarPTR;
    CmptrMemFn_void Alphabar_QG_PTR;
    z2funpointer z2limit_PTR;
    subterm_z2upper_fp z2upper_qg_PTR;
    xrapidity_y_eta_funpointer Xrpdty_LO_PTR, Xrpdty_DIP_PTR; // point to probe or target rapidity functions
    xrapidity_funpointer Xrpdty_LO_projectileY_z2min_PTR; // points to z2min consistent probe rapidity
    xrapidity_NLO_funpointer Xrpdty_NLO_PTR;
    bkkernel_funpointer K_kernel_PTR;
    trbk_rho_funpointer rho_PTR;
    INLOqg_subterm_fp ILNLOqg_subterm_PTR, ITNLOqg_subterm_PTR;
    // Cuba integrators
    // old int way: static const int vegas = 1, suave = 2, divonne = 3;
    // int cubamethod;
    string cubamethod;// Currently integrates everything with this. ->? Integrate cumbersome QG AND Bound Dipole Loop contrb. separate one?


    /*
    **  HELPERS & POINTERS
    */
    double Sr(double r, double x);
    double SrY(double r, double y);
    double SrTripole(double x01, double x_x01, double x02, double x_x02, double x21, double x_x21);
    double P(double z);
    double heaviside_theta(double x);

    // running couplings
    double Alphabar(double rsq){ return (this->*AlphabarPTR)(rsq); } // pointer shell function
    double Alphabar_QG(void *userdata){ return (this->*Alphabar_QG_PTR)(userdata); } // pointer shell function
    double alpha_bar_running_pd( double rsq );
    double alpha_bar_running_pd_sharpcutoff( double rsq );
    double alpha_bar_fixed( double rsq );

    // QG term specific couplings
    double alpha_bar_QG_running_pd( void *userdata );
    double alpha_bar_QG_fixed( void *userdata );
    double alpha_bar_QG_running_guillaume( void *userdata );
    double alpha_bar_QG_running_smallest( void *userdata );

    // z2 lower bounds
    // pointer, only one since lower bound is set consistently in all terms
    double z2lower_bound( double x, double qsq ) { return (this->*z2limit_PTR)(x,qsq); } // pointer shell function
    double z2bound_simple( double x , double qsq ) { return x/icX0 ; }
    double z2bound_improved( double x , double qsq ) { return (x/icX0)*(icQ0sqr/qsq)  ; }

    // Evolution variables / rapidities
    double Xrpdty_LO( double x, double qsq, double rsq) { return (this->*Xrpdty_LO_PTR)(x,qsq,rsq);}
    double Xrpdty_DIP( double x, double qsq, double rsq) { return (this->*Xrpdty_DIP_PTR)(x,qsq,rsq);}
    double Xrpdty_LO_simple( double x , double qsq ) { return x ; } // X = x0*z2min; z2min = xbj/x0
    double Xrpdty_LO_improved( double x , double qsq ) { return x*icQ0sqr/qsq ; } // X = x0*z2min; z2min = xbj/x0*Q0²/Q²
    double Xrpdty_NLO(double Qsq, double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq ){
        return (this->*Xrpdty_NLO_PTR)(Qsq,z2,z2min,icX0,x01sq,x02sq,x21sq);
    };

    double Xrpdty_LO_projectileY(double x, double Qsq, double rsq);
    double Xrpdty_LO_projectileY_max(double x, double Qsqr, double ){return x*icQ0sqr/Qsqr; } /// x0 = 1   
    double Xrpdty_LO_targetETA(double x, double Qsq, double rsq);

    double Xrpdty_NLO_projectileY(double Qsq, double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq );
    double Xrpdty_NLO_targetETA(double Qsq, double z2, double z2min, double icX0, double x01sq, double x02sq, double x21sq );

    // Target rapidity Eta shift calculator(s)
    double cov_to_eta_x(double r, double xbj, double Qsq);
    double rho_rapidity_shift(double Qsq, double x01sq = 0, double x02sq = 0, double x21sq = 0){ return (this->*rho_PTR)(Qsq, x01sq, x02sq, x21sq); };
    double rho_no_shift_projectileY(double Qsq, double x01sq = 0, double x02sq = 0, double x21sq = 0);
    double rho_rapidity_shift_QQ0(double Qsq, double x01sq = 0, double x02sq = 0, double x21sq = 0);
    double rho_rapidity_shift_RQ0(double Qsq, double x01sq = 0, double x02sq = 0, double x21sq = 0);
    double rho_rapidity_shift_XR(double Qsq, double x01sq, double x02sq, double x21sq);
    double rho_rapidity_shift_MAX_XYR(double Qsq, double x01sq, double x02sq, double x21sq);
    double x_eta_delta_ij_r(double x_ij_sq, double rsq);


    // BK kernels
    double K_kernel( double rsq, double x02sq, double x21sq ) { return (this->*K_kernel_PTR)(rsq,x02sq,x21sq); } // pointer shell function
    double K_resum(double rsq, double x02sq, double x21sq);
    double K_lobk(double rsq, double x02sq, double x21sq);

    // Sub scheme subtraction term pointer and its targets
    double ILNLOqg_subterm(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
        return (this->*ILNLOqg_subterm_PTR)(Q,x,z1,z2,x01sq,x02sq,x21sq);
    }
    double ILNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ILNLOqg_subterm_trbk_edmond(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    
    double ITNLOqg_subterm(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq){
        return (this->*ITNLOqg_subterm_PTR)(Q,x,z1,z2,x01sq,x02sq,x21sq);
    }
    double ITNLOqg_subterm_lobk_z2tozero(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_lobk_explicit(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_resumbk(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_kcbk_beuf(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);
    double ITNLOqg_subterm_trbk_edmond(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq);

    // Sub scheme subtraction term z2 upper bound pointer and targets
    double z2upper_qg_subterm(double z1){ return (this->*z2upper_qg_PTR)(z1); }   // pointer method
    double z2upper_unity(double z1){ return 1.0; }
    double z2upper_min_z_qqbar(double z1){
        // set z2 integration upper bound to min(z1, 1 - z1)
        double min_momentum;
        if (z1 < 0.5){
            min_momentum = z1;
        } else {
            min_momentum = 1 - z1;
        }
        return min_momentum;
    }


    // cross section integrators
    // Longitudinal
    double ILLO(double Q, double z1, double x01sq) ;
    double ILdip(double Q, double z1, double x01sq) ;
    double Bessel0Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ILNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ILNLOqgRisto(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ILNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;

    double LLOp(double Q, double x) ;
    double LLOpMass(double Q, double x, bool charm) ;
    double LLOp_massive(double Q, double x, double mf) ;
    double LNLOdip(double Q, double x) ;
    double LNLOdip_z2(double Q, double x) ;
    double LNLOdip_massive_LiLogConst(double Q, double x, double mf) ;
    double LNLOdip_massive_Iab(double Q, double x, double mf) ;
    double LNLOdip_massive_Icd(double Q, double x, double mf) ;
    double LNLOqgunsub(double Q, double x) ;
    double LNLOsigma3(double Q, double x) ;
    double LNLOqgsub(double Q, double x) ;
    double LNLOqgunsubRisto(double Q, double x) ;
    double LNLOqgsubRisto(double Q, double x) ;
    double LNLOqgunsub_massive(double Q, double x, double mf) ;
    double LNLOqgunsub_massive_I1(double Q, double x, double mf) ;
    double LNLOqgunsub_massive_I2(double Q, double x, double mf) ;
    double LNLOqgunsub_massive_I3(double Q, double x, double mf) ;


    // Transverse
    double ITLO(double Q, double z1, double x01sq) ;
    double ITdip(double Q, double z1, double x01sq) ;
    double Bessel1Tripole(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ITNLOqg(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ITNLOqgRisto(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;
    double ITNLOsigma3(double Q, double x, double z1, double z2, double x01sq, double x02sq, double x21sq) ;


    double TLOp(double Q, double x) ;
    double TLOpMass(double Q, double x, bool charm) ;
    double TLOp_massive(double Q, double x, double mf) ;
    double TNLOdip(double Q, double x) ;
    double TNLOdip_z2(double Q, double x) ;
    double TNLOdip_massive_I1(double Q, double x, double mf) ;
    double TNLOdip_massive_I2(double Q, double x, double mf) ;
    double TNLOdip_massive_I3(double Q, double x, double mf) ;
    double TNLOqgunsub(double Q, double x) ;
    double TNLOsigma3(double Q, double x) ;
    double TNLOqgsub(double Q, double x) ;
    double TNLOqgunsubRisto(double Q, double x) ;
    double TNLOqgsubRisto(double Q, double x) ;
    double TNLOqgunsub_massive_I1(double Q, double x, double mf) ;
    double TNLOqgunsub_massive_I2(double Q, double x, double mf) ;
    double TNLOqgunsub_massive_I3(double Q, double x, double mf) ;
};

double sumef_from_mass(double mf);

void Cuba(string method, int ndim, integrand_t integrand,void *userdata, double *integral, double *error, double *prob);

int integrand_ILLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip_z2(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgsub(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsubRisto(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgsubRisto(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;

int integrand_ILdip_massive_LiLogConst(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip_massive_Iab(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILdip_massive_Icd(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub_massive(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub_massive_I1(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub_massive_I2(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub_massive_I3(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub_massive_I2_fast(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ILqgunsub_massive_I3_fast(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;


int integrand_ITLOp(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITLOpMass(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITdip(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITdip_z2(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITqgunsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITsigma3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgsub(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITqgunsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;
int integrand_ITqgsubRisto(const int *ndim, const double x[], const int *ncomp, double *f, void *userdata) ;

int integrand_ITdip_massive_I1(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITdip_massive_I2(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITdip_massive_I3(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgunsub_massive_I1(const int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgunsub_massive_I2(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgunsub_massive_I3(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgunsub_massive_I2_fast(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;
int integrand_ITqgunsub_massive_I3_fast(int *ndim, const double x[], const int *ncomp,double *f, void *userdata) ;

string PrintVector(vector<double> v);

#endif // _NLODIS_SIGMAR_H
