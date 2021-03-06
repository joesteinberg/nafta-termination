#ifndef __MAIN_C__
#define __MAIN_C__

#include "globals.h"
#include "calibrate.h"
#include "eqm.h"

uint parse_args(int argc, char **argv)
{
  int opt = 0;
  uint cnt=0;

  fprintf(logfile,KBLU "Model scenario:" RESET);
  while((opt = getopt(argc, argv, "kltgjspefibanuqcdovwxr")) != -1){
    switch(opt){
    case 'k':
      k_adj_cost=0;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " No capital adjustment costs" RESET);
      break;
    case 'l':
      l_adj_cost=0;
      sens=1;
      cnt = cnt+1;
      if(cnt>1)
	cnt=1;
      fprintf(logfile,KBLU " No labor adjustment costs" RESET);
      break;
    case 't':
      f_adj_cost=0;
      m_adj_cost=0;
      sens=1;
      cnt = cnt+1;
      if(cnt>1)
	cnt=1;
      fprintf(logfile,KBLU " No trade adjustment costs" RESET);
      break;
    case 'g':
      eqkappa=1;
      sens=1;
      cnt = cnt+1;
      fprintf(logfile,KBLU " Static export participation" RESET);
      break;
    case 'j':
      nokappa=1;
      sens=1;
      cnt = cnt+1;
      if(cnt>1)
	cnt=1;
      fprintf(logfile,KBLU " No extensive margin" RESET);
      break;
    case 's':
      sym_te_flag=1;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Symmetric Armington elasticities" RESET);
      break;
    case 'p':
      ltp_flag=1;
      sens=1;
      read_seed= 1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " LTP Armington elasticities" RESET);
      break;
    case 'e':
      cobb_douglas_flag=1;
      sens=1;
      cnt = cnt+1;
      fprintf(logfile,KBLU " Cobb-Douglas production" RESET);
      break;
    case 'f':
      cobb_douglas_flag2=1;
      sens=1;
      cnt = cnt+1;
      fprintf(logfile,KBLU " Cobb-Douglas consumption" RESET);
      break;
    case 'i':
      noio_flag=1;
      m_adj_cost=0;
      cnt=cnt+1;
      fprintf(logfile,KBLU " No intermediate inputs" RESET);
      break;
    case 'b':
      iceberg_flag=1;
      sens=1;
      //cnt=cnt+1;
      fprintf(logfile,KBLU " Iceberg trade costs, not tariffs" RESET);
      break;
    case 'n':
      fix_tb_flag=1;
      sens=1;
      //cnt = cnt+1;
      fprintf(logfile,KBLU " Fixed trade balance" RESET);
      break;
    case 'a':
      fix_tb_flag2=1;
      sens=1;
      //cnt = cnt+1;
      fprintf(logfile,KBLU " Constant trade balance" RESET);
      break;
    case 'u':
      sens=1;
      us_ht_flag=1;
      cnt = cnt+1;
      fprintf(logfile,KBLU " Higher US tariffs" RESET);
      break;
    case 'q':
      sens=1;
      ucta_flag=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " US-Canada FTA" RESET);
      break;
    case 'c':
      sens=1;
      camta_flag=1;
      cnt = cnt+1;
      fprintf(logfile,KBLU " Canada-Mexico FTA" RESET);
      break;
    case 'd':
      sens=1;
      dom_con_flag=1;
      cnt = cnt+1;
      fprintf(logfile,KBLU "USMCA" RESET);
      break;
    case 'r':
      sens=1;
      dom_con_flag=2;
      cnt = cnt+1;
      fprintf(logfile,KBLU "Stricter domestic content requirements" RESET);
      break;
    case 'o':
      old_mfn_flag = 1;
      sens=1;
      //cnt = cnt+1;
      fprintf(logfile,KBLU " Pre-NAFTA MFN tariffs" RESET);
      break;
    case 'v':
      old_iomat_flag=1;
      sens=1;
      //cnt = cnt+1;
      fprintf(logfile,KBLU " Pre-NAFTA IO matrix" RESET);
      break;
    case 'w':
      fix_k_flag = 1;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " Fixed capital stocks" RESET);
      break;
    case 'x':
      no_k_flag = 1;
      sens=1;
      cnt=cnt+1;
      fprintf(logfile,KBLU " No capital" RESET);
      break;
    case '?':
      fprintf(logfile,"\nIncorrect command line option: %c. Possible options:\n",opt);
      fprintf(logfile,"\t-k: no capital adjustment costs\n");
      fprintf(logfile,"\t-l: no labor adjustment costs\n");
      fprintf(logfile,"\t-t: no trade adjustment costs\n");
      fprintf(logfile,"\t-g: static export participation (kappa0 = kappa1)\n");
      fprintf(logfile,"\t-j: no extensive margin (kappa0 = kappa1 = 0)\n");
      fprintf(logfile,"\t-s: symmetric trade elasticities\n");
      fprintf(logfile,"\t-p: LTP trade elasticities\n");
      fprintf(logfile,"\t-e: cobb-douglas production\n");
      fprintf(logfile,"\t-f: cobb-douglas consumption\n");
      fprintf(logfile,"\t-i: no intermediate inputs\n");
      fprintf(logfile,"\t-b: iceberg costs, not tariffs\n");
      fprintf(logfile,"\t-u: higher US tariffs\n");
      fprintf(logfile,"\t-q: USA-CAN FTA\n");
      fprintf(logfile,"\t-c: CAN-MEX FTA\n");
      fprintf(logfile,"\t-d: USMCA\n");
      fprintf(logfile,"\t-r: Stricter domestic content requirements in all sectors\n");
      fprintf(logfile,"\t-o: Pre-NAFTA MFN tariffs\n");
      fprintf(logfile,"\t-v: Pre-NAFTA IO matrix\n");
      fprintf(logfile,"\t-w: Fixed capital stocks\n");
      fprintf(logfile,"\t-x: No capital\n");
      return 1;
      break;
    }
  }

  if(cnt>1)
    {
      fprintf(logfile,KRED "\nOnly one sensitivity analysis allowed at a time!\n" RESET);
      return 1;
    }
  else
    {
      if(cnt==0 && sens==0)
	{
	  fprintf(logfile,KBLU " Baseline analysis" RESET);
	  write_seed = 1;
	}
      fprintf(logfile,"\n");
  
      return 0;
    }
}

int quant_exercise()
{
  // -----------------------------------------------------------------------------------------------------------
  // set up variable and parameter structures
  uint it;
  if(calibrate())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  
  for(it=0; it<NTH; it++)
    {
      init_vars(&(eee0[it]));
      init_vars(&(eee1[it]));
    }

  // -------------------------------------------------------------------------------------------------------
  // no repeal

  fprintf(logfile, KNRM "----------------------------------------------------------------------\n" RESET);

  scenario = 0;
  set_neqm();

  fprintf(logfile,KBLU "\nSolving for no-repeal equilibrium...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee0[0]), &(ppp0[0]));

  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_usa",0);
  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_can",1);
  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_mex",2);
  write_eqm_vars(&(eee0[0]),&(ppp0[0]),"vars0_row",3);

  // -------------------------------------------------------------------------------------------------------
  // with repeal
  scenario = 1;
  set_neqm();
  for(it=0; it<NTH; it++)
    {
      set_tariffs(&(ppp1[it]),scenario);
    }

  fprintf(logfile,KBLU "\nSolving for equilibrium with repeal...\n" RESET);
  if(solve_eqm())
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    }
  calc_welfare(&(eee1[0]), &(ppp1[0]));

  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_usa",0);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_can",1);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_mex",2);
  write_eqm_vars(&(eee1[0]),&(ppp1[0]),"vars1_row",3);

  // -------------------------------------------------------------------------------------------------------
  // 1-period average SR trade elasticity
  double tmp = (eee1[0].te_t[TNAFTA][0] + eee1[0].te_t[TNAFTA][1] + eee1[0].te_t[TNAFTA][2])/3.0;
  fprintf(logfile,KBLU "\nAverage SR trade elasticity: %0.4f\n" RESET,tmp);
  
  return 0;
}

int main(int argc, char * argv[])
{
  noio_flag=0;
  sens=0;
  dom_con_flag=0;
  par = 0;
  iceberg=1;
  solver_verbose=1;
  cobb_douglas_flag=0;
  cobb_douglas_flag2=0;
  f_adj_cost=1;
  m_adj_cost=1;
  k_adj_cost=1;
  l_adj_cost=1;
  fixl=1;
  iceberg_flag=0;
  eval_eqm_once_flag=0;
  eval_bgp_once_flag=0;
  read_seed=0;
  write_seed=0;
  camta_flag=0;
  ucta_flag=0;
  sym_te_flag=0;
  us_ht_flag=0;
  us_ht_flag2=0;
  ltp_flag=0;
  old_mfn_flag=0;
  fix_tb_flag=0;
  fix_tb_flag2=0;
  old_iomat_flag=0;
  eqkappa=0;
  nokappa=0;
  fix_k_flag=0;
  no_k_flag=0;
  eval_calfn_once=0;
  logfile = stdout;

  fprintf(logfile, KGRN "\nThe Macroeconomic Effects of NAFTA Repeal" RESET);
  fprintf(logfile, KGRN "\nJoseph Steinberg, University of Toronto" RESET);
  fprintf(logfile, KNRM "\n" RESET);

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KBLU "\nSetting up environment...\n\n" RESET);

  if(parse_args(argc,argv))
    {
      fprintf(logfile, KRED "\nProgram failed!\n" RESET);
      return 1;
    };

  // -----------------------------------------------------------------------
  // set up parallel environment
#ifdef _OPENMP
  omp_set_num_threads(NTH);
  mkl_set_num_threads(NTH);
  uint nt = omp_get_max_threads();
  fprintf(logfile, KBLU "\n\tParallel processing using %d OMP threads\n" RESET,nt);
#pragma omp parallel num_threads(nt)
  {
    int it = omp_get_thread_num();
    fprintf(logfile,KBLU "\t\tHello from thread %d out of %d\n" RESET,
	    it, nt);
  }
  fprintf(logfile,"\n");
#endif

  quant_exercise();

  // ------------------------------------------------------------------------

  fprintf(logfile, KNRM "\n----------------------------------------------------------------------\n" RESET);
  fprintf(logfile, KGRN "\nProgram complete!\n\n" RESET);

  return 0;
}

#endif
