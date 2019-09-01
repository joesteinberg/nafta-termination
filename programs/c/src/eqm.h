#ifndef __EQM_H__
#define __EQM_H__

#include "globals.h"
#include "calibrate.h"
#include "solver.h"

extern const uint nbgp; // number of BGP variables... used to solve for BGP to create initial guess for equilibrium
uint neqm; // number of equilibrium vars/eqns... depends on situation
double bbgp[NC-1]; // BGP bond-holdings... state variable
uint scenario; // scenario... may be redundant in stochastic model

// eqm contains all the vars associated with a deterministic equilibrium, or a single history of a stochastic equilibrium
typedef struct
{
  double pb_t[NT+1];

  double b_t[NT+1][NC];
  double cc_t[NT+1][NC];
  double ii_t[NT+1][NC];
  double ll_t[NT+1][NC];
  double kk_t[NT+1][NC];
  double cpi_t[NT+1][NC];
  double pi_t[NT+1][NC];
  double w_t[NT+1][NC];
  double rk_t[NT+1][NC];
  double ngdp_t[NT+1][NC];
  double ngdp2_t[NT+1][NC];
  double rgdp_t[NT+1][NC];
  double iy_t[NT+1][NC];
  double lp_agg_t[NT+1][NC];

  double ex_t[NT+1][NC][NC];
  double im_t[NT+1][NC][NC];
  double nx_t[NT+1][NC][NC];
  double exf_t[NT+1][NC][NC];
  double imf_t[NT+1][NC][NC];
  double nxf_t[NT+1][NC][NC];
  double exm_t[NT+1][NC][NC];
  double imm_t[NT+1][NC][NC];
  double nxm_t[NT+1][NC][NC];
  double rer_t[NT+1][NC][NC];

  double rex_t[NT+1][NC][NC];
  double rim_t[NT+1][NC][NC];
  double rexf_t[NT+1][NC][NC];
  double rimf_t[NT+1][NC][NC];
  double rexm_t[NT+1][NC][NC];
  double rimm_t[NT+1][NC][NC];

  double y_t[NT+1][NC][NS];
  double MC_t[NT+1][NC][NS];
  double ws_t[NT+1][NC][NS];
  double rks_t[NT+1][NC][NS];
  double py_t[NT+1][NC][NS];
  double va_t[NT+1][NC][NS];
  double md_t[NT+1][NC][NS][NS];
  double k_t[NT+1][NC][NS];
  double kd_t[NT+1][NC][NS];
  double l_t[NT+1][NC][NS];
  double lf_t[NT+1][NC][NS];
  double is_t[NT+1][NC][NS];
  double lp_t[NT+1][NC][NS];

  double exs_t[NT+1][NC][NS][NC];
  double ims_t[NT+1][NC][NS][NC];
  double nxs_t[NT+1][NC][NS][NC];
  double rexs_t[NT+1][NC][NS][NC];
  double rims_t[NT+1][NC][NS][NC];

  double pm_t[NT+1][NC][NS];
  double m_t[NT+1][NC][NS];  
  double m2_t[NT+1][NC][NS][NC];
  double m2s_t[NT+1][NC][NS][NC];
  double pm2_t[NT+1][NC][NS][NC];

  double p_t[NT+1][NC][NS];
  double q_t[NT+1][NC][NS];  
  double q2_t[NT+1][NC][NS][NC];
  double q2s_t[NT+1][NC][NS][NC];
  double pq2_t[NT+1][NC][NS][NC];

  double c_t[NT+1][NC][NS];
  double i_t[NT+1][NC][NS];

  double welfare_t[NT+1][NC];
  double welfare2_t[NT+1][NC];
  double welfare3_t[NT+1][NC];
  double welfare4_t[NT+1][NC];
  double welfare_cost_t[NT+1][NC];

  double Q_t[NT+1][NC];

  double aes2_m_t[NT+1][NC][NS][NC];
  double aes2_f_t[NT+1][NC][NS][NC];
  double aes2_t[NT+1][NC][NS][NC];
  double aes_t[NT+1][NC][NS];
  double ae_t[NT+1][NC];

  double tes2_t[NT+1][NC][NS][NC];
  double tes_t[NT+1][NC][NS];
  double te_t[NT+1][NC];
  double realloc2_t[NT+1][NC][NS];
  double realloc_t[NT+1][NC];

  double Df_t[NT+1][NC][NS][NC];
  double Dbf_t[NT+1][NC][NS][NC];
  double Dhf_t[NT+1][NC][NS][NC];
  double Dtf_t[NT+1][NC][NS][NC];

  double Dm_t[NT+1][NC][NS][NC];
  double Dbm_t[NT+1][NC][NS][NC];
  double Dhm_t[NT+1][NC][NS][NC];
  double Dtm_t[NT+1][NC][NS][NC];

  double exrate_t[NT+1][NC][NS][NC];
  double exitrate_t[NT+1][NC][NS][NC];
  
  exporter_vars ev_t[NT+1][NC][NS][NC-1];
}eqm;

// array of NTH eqm structs for use in solving deterministic equilibrium (like no-Brexit counterfactual)
// we use the array to parallelize the process of evaluating the jacobian matrix in the solver
eqm eee0[NTH];
eqm eee1[NTH];

void set_neqm(); // sets the dimension of the equilibrium solution space
void init_vars(eqm * e); // initializes all the variables of an eqm struct to zero (or 1 where appropriate)
void copy_vars(eqm * e1, const eqm * e0); // copies all the variables from one eqm struct to another
uint stack_bgp_vars(double * myx, const eqm * e); // stacks the BGP variables (last period of an eqm struct) into an array
uint unstack_bgp_vars(eqm * e, const double * myx); // unstacks BGP vars from array to last period of an eqm struct
uint stack_eqm_vars(double * myx, const eqm * e); // stacks deterministic equilibrium vars into an array
uint unstack_eqm_vars(eqm * e, const double * myx); // unstacks deterministic equillibrium vars
uint set_initial_bpg_guess(); // constructs initial guess for a BGP... the "initial guess for the initial guess" function
uint write_eqm_vars(const eqm * e, const params * p, char * fname, uint i); // write main deterministic equilibrium vars for country i to file
void set_vars(eqm * e, const params * p, uint t, uint bgp); // sets all the variables for a given period t
uint eval_bgp_conds(const double * myx, double * myf, uint tn); // evaluates the BGP equations
uint solve_bgp(double bb[NC-1]); // solves for the balanced growth path
uint eval_eqm_conds(const double * myx, double * myf, uint tn); // evaluates the deterministic equilibrium conditions
uint solve_eqm(); // solves for the deterministic equilibrium
void calc_welfare(eqm * e, const params * p);

// inlined equilibrium equations
static inline double mkt_clear_k(const params * p, const eqm * e, uint t, uint i, uint s)
{
  //double retval = (e->k_t[t][i][s] - e->kd_t[t][i][s])/e->k_t[t][i][s];
  double retval = (e->k_t[t][i][s] - e->kd_t[t][i][s]);
  return retval;
}

static inline double mkt_clear_q2(const params * p, const eqm * e, uint t, uint i, uint s, uint j)
{
  //double retval = (e->q2_t[t][i][s][j] - e->q2s_t[t][i][s][j])/e->q2_t[t][i][s][j];
  double retval = (e->q2_t[t][i][s][j] - e->q2s_t[t][i][s][j]);
  return retval;
}

static inline double mkt_clear_m2(const params * p, const eqm * e, uint t, uint i, uint s, uint j)
{
  //double retval = (e->m2_t[t][i][s][j] - e->m2s_t[t][i][s][j])/e->m2_t[t][i][s][j];
  double retval = (e->m2_t[t][i][s][j] - e->m2s_t[t][i][s][j]);
  
  return retval;
}

static inline double mkt_clear_m(const params * p, const eqm * e, uint t, uint i, uint s)
{
  double retval = e->m_t[t][i][s];
  uint r;
  for(r=0; r<NS; r++)
    {
      retval = retval - e->md_t[t][i][r][s];
    }
  //retval = retval/e->m_t[t][i][s];
  return retval;
}

static inline double mkt_clear_q(const params * p, const eqm * e, uint t, uint i, uint s)
{
  //return (e->q_t[t][i][s] - e->c_t[t][i][s] - e->i_t[t][i][s])/e->q_t[t][i][s];
  return (e->q_t[t][i][s] - e->c_t[t][i][s] - e->i_t[t][i][s]);
}

static inline double mucs_mucr(const params * p, const eqm * e, uint t, uint i, uint s, uint r)
{
  return p->eps[i][0][s]*pow(e->c_t[t][i][s],p->rho-1.0)/e->p_t[t][i][s]
    - p->eps[i][0][r]*pow(e->c_t[t][i][r],p->rho-1.0)/e->p_t[t][i][r];
}

static inline double muc_mul(const params * p, const eqm * e, uint t, uint i)
{
  if(fixl==1)
    {
      //return (e->ll_t[t][i] - p->lbar[i]/3.0)/(p->lbar[i]/3.0);
      return (e->ll_t[t][i] - p->lbar[i]/3.0);
    }
  else
    {
      return 1.0 * (
		    muc
		    (
		     e->c_t[t][i],
		     e->ll_t[t][i],
		     p->lbar[i],
		     p->eps[i][0],
		     p->rho,
		     p->phi[i],
		     p->psi,
		     2)/e->p_t[t][i][2] - 
		    mul
		    (
		     e->c_t[t][i],
		     e->ll_t[t][i],
		     p->lbar[i],
		     p->eps[i][0],
		     p->rho,
		     p->phi[i],
		     p->psi)/e->w_t[t][i] );
    }
}

static inline double euler(const params * p, const eqm * e, uint t, uint i)
{
  if(fix_tb_flag && scenario>0 && i<(NC-1))
    {
      return e->b_t[t+1][i] - eee0[0].b_t[t+1][i];
    }
  else if(fix_tb_flag2 && scenario>0 && i<(NC-1))
    {
      return e->b_t[t+1][i] - eee0[0].b_t[t+1][i];
    }
  else
    {
      return 10000.0 * (e->pb_t[t] * muc(
					 e->c_t[t][i],
					 e->ll_t[t][i],
					 p->lbar[i],
					 p->eps[i][0],
					 p->rho,
					 p->phi[i],
					 p->psi,
					 2)
			- p->beta[i] * e->cpi_t[t+1][0] * (e->p_t[t][i][2]/e->p_t[t+1][i][2])
			* muc(
			      e->c_t[t+1][i],
			      e->ll_t[t+1][i],
			      p->lbar[i],
			      p->eps[i][0],
			      p->rho,
			      p->phi[i],
			      p->psi,
			      2));
    }
}

static inline double bop(const params * p, const eqm * e, uint t, uint i)
{
  double net_iceberg = 0.0;
  if(iceberg==1 && 0)
    {
      uint j;
      for(j=0; j<NC; j++)
	{
	  if(j!=i)
	    {
	      uint s;
	      for(s=0; s<NS; s++)
		{
		  net_iceberg = net_iceberg + 
		    + e->pm2_t[t][j][s][j] * p->ntb_m_ts[t][j][s][i]*e->m2_t[t][j][s][i] 
		    + e->pq2_t[t][j][s][i] * p->ntb_f_ts[t][j][s][i]*e->q2_t[t][j][s][i]
		    - e->pm2_t[t][i][s][j] * p->ntb_m_ts[t][i][s][j]*e->m2_t[t][i][s][j] 
		    - e->pq2_t[t][i][s][j] * p->ntb_f_ts[t][i][s][j]*e->q2_t[t][i][s][j];
		}
	    }
	}
    }

  if(t<NT)
    {
      return SUM(e->nx_t[t][i],NC) + net_iceberg + e->b_t[t][i]*e->cpi_t[t][0] - e->pb_t[t]*e->b_t[t+1][i];
    }
  else
    {
      return SUM(e->nx_t[t][i],NC) + net_iceberg + e->b_t[t][i]*e->cpi_t[t][0] - e->pb_t[t]*e->b_t[t][i];
    }
}

static inline double price_norm(const eqm * e, uint t)
{
  return e->cpi_t[t][0] - 1.0;
}

static inline double foc_m2(const params * p, const eqm * e, uint t, uint i, uint s, uint j)
{
  if(noio_flag)
    {
      return e->m2_t[t][i][s][j];
    }
  else
    {
      double tmp = e->pm_t[t][i][s] * p->mu[i][s][j] * pow(p->M[i][s],p->zeta[i][s])
	* pow(e->m_t[t][i][s]/e->m2_t[t][i][s][j],1.0 - p->zeta[i][s]);

      tmp = tmp - (1.0+p->tau_m_ts[t][i][s][j]) * e->pm2_t[t][i][s][j];

      //if(t<NT && m_adj_cost)
      if(t<NT)
	{
	  if(j != i)
	    {
	      if(t>0)
		{
		  tmp = tmp - e->w_t[t][i] * p->etaM * (2.0 * e->m2_t[t][i][s][j]/e->m2_t[t-1][i][s][j] - 2.0);
		}
	      else
		{
		  tmp = tmp - e->w_t[t][i] * p->etaM * (2.0 * e->m2_t[t][i][s][j]/p->m02[i][s][j] - 2.0);
		}

	      if(t<NT-1)
		tmp = tmp - e->Q_t[t][i] * e->w_t[t+1][i] * p->etaM * (1.0 - e->m2_t[t+1][i][s][j]*e->m2_t[t+1][i][s][j]/e->m2_t[t][i][s][j]/e->m2_t[t][i][s][j]);
      
	    }
	}

      return tmp;
    }
}

static inline double foc_q2(const params * p, const eqm * e, uint t, uint i, uint s, uint j)
{
  double tmp = e->p_t[t][i][s] * p->theta[i][s][j] * pow(p->H[i][s],p->sig[i][s])
    * pow(e->q_t[t][i][s]/e->q2_t[t][i][s][j],1.0 - p->sig[i][s]);

  tmp = tmp - (1.0+p->tau_f_ts[t][i][s][j]) * e->pq2_t[t][i][s][j];

  //if(t<NT && f_adj_cost)
  if(t<NT)
    {
      if(j != i)
	{
	  if(t>0)
	    {
	      tmp = tmp - e->w_t[t][i] * p->etaF * (2.0 * e->q2_t[t][i][s][j]/e->q2_t[t-1][i][s][j] - 2.0);
	    }
	  else
	    {
	      tmp = tmp - e->w_t[t][i] * p->etaF * (2.0 * e->q2_t[t][i][s][j]/p->q02[i][s][j] - 2.0);
	    }

	  if(t<NT-1)
	    {
	      tmp = tmp - e->Q_t[t][i] * e->w_t[t+1][i] * p->etaF * (1.0 - e->q2_t[t+1][i][s][j]*e->q2_t[t+1][i][s][j]/e->q2_t[t][i][s][j]/e->q2_t[t][i][s][j]);
	    }
	}
    }

  return tmp;
}

static inline double noarb(const params * p, const eqm * e, uint t, uint  i, uint s)
{
  double retval = 0.0;
  
  if(!fix_k_flag)
    {
      retval = e->rks_t[t][i][s] - e->rk_t[t][i];
    
      //  if(k_adj_cost && t<NT)
      //if(t<NT)
      {
	retval -= e->w_t[t][i] * p->etaK * (2.0 * e->k_t[t][i][s]/e->k_t[t-1][i][s] - 2.0);
	if(t<(NT))
	  {
	    retval -= e->Q_t[t][i]*e->w_t[t+1][i]*p->etaK*
	      (1.0-e->k_t[t+1][i][s]*e->k_t[t+1][i][s]/e->k_t[t][i][s]/e->k_t[t][i][s]);
	  }
      }
    }
  else
    {
      retval = e->k_t[t][i][s] - p->k0[i][s];
    }
  
  return retval; 
}

static inline double noarb_l(const params * p, const eqm * e, uint t, uint  i, uint s)
{
  double retval = e->ws_t[t][i][s] - e->w_t[t][i];

  if(t==0)
    {
      retval -= e->w_t[t][i] * p->etaL * (2.0 * e->l_t[t][i][s]/p->l0[i][s] - 2.0);
    }
  else if(t<NT)
    {
      retval -= e->w_t[t][i] * p->etaL * (2.0 * e->l_t[t][i][s]/e->l_t[t-1][i][s] - 2.0);
    }
  if(t<(NT))
    {
      retval -= e->Q_t[t][i]*e->w_t[t+1][i]*p->etaL*
	(1.0-e->l_t[t+1][i][s]*e->l_t[t+1][i][s]/e->l_t[t][i][s]/e->l_t[t][i][s]);
    }
  
  return retval; 
}

int bgp_func_f(const gsl_vector * x, void * data, gsl_vector * f);
int bgp_func_df(const gsl_vector * x, void * data, gsl_matrix * J);
int bgp_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);
int eqm_func_f(const gsl_vector * x, void * data, gsl_vector * f);
int eqm_func_df(const gsl_vector * x, void * data, gsl_matrix * J);
int eqm_func_fdf(const gsl_vector * x, void * data, gsl_vector * f, gsl_matrix * J);

#endif
