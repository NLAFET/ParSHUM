#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>

#include "TP_verbose.h"
#include "TP_enum.h"
#include "TP_auxiliary.h"

TP_verbose_parms
TP_verbose_default_parms()
{
  TP_verbose_parms self = calloc(1, sizeof(*self));
  char *output_dir = malloc( (PATH_LENGTH + 1) * sizeof(*output_dir)); 

  output_dir = getcwd(output_dir, (PATH_LENGTH + 1) * sizeof(*output_dir));

  self->prog_name       = "TP_solver";
  self->output_dir      = output_dir;
  self->user_out_dir    = 0;
  self->outfiles_prefix = NULL;

  self->out_file  = stdout;
  self->user_out_file = 0;
  self->verbosity = 1;

  return self;
}

TP_verbose 
TP_verbose_create_V0(TP_exe_parms exe_parms) 
{
  TP_verbose self  = calloc((size_t) 1, sizeof(*self));
  
  self->parms     = TP_verbose_default_parms();
  self->exe_parms = exe_parms;
  self->reason    = TP_reason_unknown;
  
  return self;
}

TP_verbose_per_step
TP_verbose_step_start_V0(TP_verbose self)
{
  TP_verbose_per_step step = calloc((size_t) 1, sizeof(*step));
  
  if (!self->stats_first) {
    self->stats_first = self->stats_last = step;
  } else {
    self->stats_last->next = step;
    self->stats_last = step;
  }
  
  self->nb_steps++;
  return step;
}

void
TP_verbose_update_Luby_step_V0(TP_verbose_per_step step, TP_Luby_step Luby_step)
{
  if (!step->nb_Luby_steps) 
    {
      step->Luby_step_first = step->Luby_step_last = Luby_step;
    }
  else
    {
      step->Luby_step_last->next = Luby_step;
      step->Luby_step_last = Luby_step;
    }
  step->nb_Luby_steps++;
}

void
TP_verbose_print_steps(TP_verbose_per_step self, TP_verbose_parms parms)
{
  FILE *file      = parms->out_file;
  char *prog_name = parms->prog_name;

  fprintf(file,"[%s]\n", prog_name);
  fprintf(file,"[%s]----------------------------------------------------------------------------------------------------------------------------------------\n", prog_name);
  fprintf(file,"[%s] ( nb_pivots )\t|| Pivot search\t| extracting\t| merging pivot\t|| Update Schur\t|   LD\t\t|   U\t\t|    S\t\t||   \n", prog_name);
  fprintf(file,"[%s]----------------------------------------------------------------------------------------------------------------------------------------\n", prog_name);
  while(self)
    {
      fprintf(file,"[%s] (%4d)  \t\t|| %12f\t| %12f\t| %12f\t|| %12f\t| %12f\t| %12f\t| %12f\t||\n", prog_name, self->nb_pivots,
	      self->timing_pivot_search, self->timing_extracting_candidates, self->timing_merging_pivots, 
	      self->timing_apply_perms, self->timing_update_LD, self->timing_update_U, self->timing_update_S);
      self = self->next;
    }
  fprintf(file,"[%s]----------------------------------------------------------------------------------------------------------------------------------------\n", prog_name);
}


void
TP_verbose_print_parms_raw(TP_exe_parms exe_parms, TP_parm_type type, FILE *file)
{
  fprintf(file,"###################PARAMETERS########################\n");
  fprintf(file,"#matrix\t%s\n", exe_parms->matrix_file);
  if (type != TP_value_tol)
    fprintf(file,"#valut_tol\t%f\n", exe_parms->value_tol);
  if (type != TP_marko_tol)
    fprintf(file,"#marko_tol\t%f\n", exe_parms->marko_tol);
  fprintf(file,"#extra_space\t%f\n", exe_parms->extra_space);
  fprintf(file,"#extra_space_inbetween\t%f\n", exe_parms->extra_space_inbetween);
  if (type != TP_schur_density)
    fprintf(file,"#density_tolerance\t%f\n", exe_parms->density_tolerance);
  if (type != TP_min_pivots)
    fprintf(file,"#min_pivot_per_steps\t%d\n", exe_parms->min_pivot_per_steps);
  if (type != TP_nb_threads)
    fprintf(file,"#nb_threads\t%d\n", exe_parms->nb_threads);
  if (type != TP_nb_candidates)
    fprintf(file,"#nb_candidates_per_block\t%d\n", exe_parms->nb_candidates_per_block);
  fprintf(file,"#nb_previous_pivots\t%d\n", exe_parms->nb_previous_pivots);
  fprintf(file,"#\n");
}

void
TP_verbose_print_group_run(TP_verbose verbose, TP_parm_type type,
			   void *val, int current_run, FILE *file)
{
  if( !current_run ) {
    switch (type) {
    case (TP_value_tol) :
      fprintf(file, "\"Value tolerance\"\t");
      break;
    case (TP_marko_tol) :
      fprintf(file, "\"Markowitz factor\"\t");
      break;
    case (TP_schur_density) :
      fprintf(file, "\"Schur density\"\t");
      break;
    case (TP_nb_candidates) :
      fprintf(file, "\"#candidates\"\t");
      break;
    case (TP_min_pivots) :
      fprintf(file, "\" #min pivots\"\t");
      break;
    case (TP_nb_threads) :
      fprintf(file, "\" #threads\"\t");
      break;
    default :
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "unrecognized type of exe_parms");
    }
    
    fprintf(file, " \"timing total\" \"timing sparse\"  \"timing dense\" \"timing convert\"");
    fprintf(file, " \"timing pivot search\" \"timing creating init sets\"  \"timing merging pivots\"");
    fprintf(file, " \"timing update\" \"timing update LD\"  \"timing update U\" \"timing update S\"");
    fprintf(file, " \"timing solve\" \"timing solve L\"  \"timing solve dense\" \"timing solve U\"");
    fprintf(file, " \"forward error\" \"backward error\"  \"schur density\" \"sparse pivots S\"");
    fprintf(file, " \"dense pivots\" \"nnz\"  \"nnz L\" \"nnz U \" #steps\n");
  }
  
  switch (type) {
  case (TP_value_tol) :
  case (TP_marko_tol) :
  case (TP_schur_density) :
    fprintf(file, "%f\t", *((double *) val));
    break;
  case (TP_nb_candidates) :
  case (TP_min_pivots) :
  case (TP_nb_threads) :
    fprintf(file, "%d\t", *((int *) val));
    break;
  default :
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "unrecognized type of exe_parms");
  }
  
  fprintf(file, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%e\t%f\t%d\t%d\t%ld\t%ld\t%ld\t%d\n",
	  verbose->timing_facto,
	  verbose->timing_facto_sparse,
	  verbose->timing_facto_dense,
	  verbose->timing_convert_schur,
	  verbose->timing_total_pivot_search,
	  verbose->timing_total_extracting_candidates, 
	  verbose->timing_total_merging_pivots,
	  verbose->timing_total_update, 
	  verbose->timing_total_update_LD, 
	  verbose->timing_total_update_U, 
	  verbose->timing_total_update_S,
	  verbose->timing_solve, 
	  verbose->timing_solve_L, 
	  verbose->timing_solve_dense, 
	  verbose->timing_solve_U,
	  verbose->forward_error, 
	  verbose->backward_error, 
	  verbose->schur_density, 
	  verbose->sparse_pivots,
	  verbose->dense_pivots, 
	  verbose->nnz_final, 
	  verbose->nnz_L, 
	  verbose->nnz_U, 
	  verbose->nb_steps);
}

void
TP_verbose_print_raw(TP_verbose verbose)
{
  TP_verbose_per_step step = verbose->stats_first;
  TP_exe_parms exe_parms = verbose->exe_parms;
  char filename[2048];
  FILE *file;
  double pivots_avg = 0.0, pivots_stddev = 0.0;

  snprintf(filename, 2048, "%s/data/%s_raw.dat", verbose->parms->output_dir, verbose->parms->outfiles_prefix);
  file = fopen(filename, "w+");

  while(step)
    {
      pivots_avg += step->nb_pivots;
      verbose->timing_total_pivot_search          += step->timing_pivot_search;
      verbose->timing_total_extracting_candidates += step->timing_extracting_candidates;
      verbose->timing_total_merging_pivots        += step->timing_merging_pivots;
      verbose->timing_total_update_LD             += step->timing_update_LD;
      verbose->timing_total_update_U              += step->timing_update_U;
      verbose->timing_total_update_S              += step->timing_update_S;
      verbose->timing_total_update += step->timing_update_LD + step->timing_update_U + step->timing_update_S;
      step = step->next;
    }
  step = verbose->stats_first;
  pivots_avg /= verbose->nb_steps;

  while(step)
    {
      pivots_stddev = (step->nb_pivots - pivots_avg) * (step->nb_pivots - pivots_avg);
      step = step->next;
    }
  step = verbose->stats_first;
  pivots_stddev /= verbose->nb_steps;
  pivots_stddev = sqrt(pivots_stddev);
  
  TP_verbose_print_parms_raw(exe_parms, TP_parm_none, file);
  fprintf(file,"####################OUTPUT##########################\n");
  fprintf(file,"#        ###########TIMING##################\n");
  fprintf(file,"#timing_facto\t%f\n", verbose->timing_facto);
  fprintf(file,"#timing_facto_sparse\t%f\n", verbose->timing_facto_sparse);
  fprintf(file,"#timing_facto_dense\t%f\n", verbose->timing_facto_dense);
  fprintf(file,"#timing_total_pivot_search\t%f\n", verbose->timing_total_pivot_search);
  fprintf(file,"#timing_total_merging_pivots\t%f\n", verbose->timing_total_merging_pivots);
  fprintf(file,"#timing_total_extracting_candidates\t%f\n", verbose->timing_total_extracting_candidates);
  fprintf(file,"#timing_total_update\t%f\n", verbose->timing_total_update);
  fprintf(file,"#timing_total_update_LD\t%f\n", verbose->timing_total_update_LD);
  fprintf(file,"#timing_total_update_U\t%f\n",  verbose->timing_total_update_U);
  fprintf(file,"#timing_total_update_S\t%f\n",  verbose->timing_total_update_S);
  fprintf(file,"#timing_convert_schur\t%f\n", verbose->timing_convert_schur);
  fprintf(file,"#timing_solve\t%f\n", verbose->timing_solve);
  fprintf(file,"#timing_solve_L\t%f\n", verbose->timing_solve_L);
  fprintf(file,"#timing_solve_dense\t%f\n", verbose->timing_solve_dense);
  fprintf(file,"#timing_solve_U\t%f\n", verbose->timing_solve_U);
  fprintf(file,"#        ###########SIZES###################\n");
  fprintf(file,"#n\t%d\n", verbose->n);
  fprintf(file,"#m\t%d\n", verbose->m);
  fprintf(file,"#nnz_input\t%ld\n", verbose->nnz_input);
  fprintf(file,"#nnz_final\t%ld\n", verbose->nnz_final);
  fprintf(file,"#nnz_L\t%ld\n", verbose->nnz_L);
  fprintf(file,"#nnz_U\t%ld\n", verbose->nnz_U);
  fprintf(file,"#nnz_S_dense\t%ld\n", verbose->nnz_S_dense);
  fprintf(file,"#        ###########OTHERS##################\n");  
  if (verbose->Luby) 
    fprintf(file,"#Luby_algorithm\n");
  fprintf(file,"#schur_density\t%f\n", verbose->schur_density);
  fprintf(file,"#pivots_avg\t%f\n", pivots_avg);
  fprintf(file,"#pivots_stddev\t%f\n", pivots_stddev);
  fprintf(file,"#forward_error\t%e\n", verbose->forward_error);
  fprintf(file,"#backward_error\t%e\n", verbose->backward_error);
  fprintf(file,"#sparse_pivots\t%d\n", verbose->sparse_pivots);
  fprintf(file,"#dense_pivots\t%d\n", verbose->dense_pivots);
  fprintf(file,"#nb_steps\t%d\n", verbose->nb_steps);

  if (verbose->reason & TP_reason_unknown) 
    fprintf(file,"#unknown_reason_for_switch\n");
  else if (verbose->reason & TP_reason_density) 
    fprintf(file,"#density_failed\n");
  else if (verbose->reason & TP_reason_no_pivots)  
    fprintf(file,"#min_pivots_failed\n");
  else if (verbose->reason & TP_reason_because) 
    fprintf(file,"#switched_to_dense_just_because\n");

  if (verbose->reason & TP_reason_dense_too_large)
    fprintf(file,"#facto_unfinished(too_large_dense_part)\n");

  while(step)
    {
      fprintf(file,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", step->nb_pivots,
	      step->timing_pivot_search, step->timing_extracting_candidates, step->timing_merging_pivots, 
	      step->timing_apply_perms, step->timing_update_LD, step->timing_update_U, step->timing_update_S);
      step = step->next;
    }
  fclose(file);

  TP_verbose_print_Luby(verbose);
}

void
TP_verbose_print_Luby(TP_verbose self)
{
  TP_verbose_per_step step = self->stats_first;
  char filename[2048];
  FILE *file;

  snprintf(filename, 2048, "%s/data/%s_Luby.dat", self->parms->output_dir, self->parms->outfiles_prefix);
  file = fopen(filename, "w+");

  fprintf(file,"===================================================================================================================================================<\n");
  fprintf(file,"| #candidates\t| #pivots\t| candidates\t| scores\t|computing max\t| first passs\t|  second pass\t|  discarding\t| managing pivots |\n");
  fprintf(file,"===================================================================================================================================================\n");
  while(step)
    {
      TP_Luby_step init_phase = step->Luby_init_phase;
      TP_Luby_step Luby_step = step->Luby_step_first;

      /* this test is needed because of singletons */
      if ( init_phase ) {
	fprintf(file,"| %12d\t| \t\t| %12f\t| %12f\t| 0.000000000\t| 0.000000000\t| 0.000000000\t| 0.000000000\t| %12f\t|\n", 
		init_phase->nb_candidates, init_phase->timing_first_pass,
		init_phase->timing_second_pass, init_phase->timing_discarding);
	while(Luby_step) {
	  fprintf(file,"| %12d\t| %12d\t| 0.000000000\t| 0.000000000\t| %12f\t| %12f\t| %12f\t| %12f\t| 0.000000000\t|\n",
		  Luby_step->nb_candidates, Luby_step->nb_pivots,
		  Luby_step->timing_max, Luby_step->timing_first_pass,
		  Luby_step->timing_second_pass, Luby_step->timing_discarding);
	  Luby_step = Luby_step->next;
	}
	fprintf(file,"===================================================================================================================================================\n");
      }
      step = step->next;
    }
  fclose(file);
}

void
TP_verbose_print_V0(TP_verbose self)
{
  FILE *file             = self->parms->out_file;
  char *prog_name        = self->parms->prog_name;
  TP_exe_parms exe_parms = self->exe_parms;
  
  fprintf(file,"[%s] \n", prog_name);
  fprintf(file,"[%s] Parameters:\n", prog_name);
  if (exe_parms->matrix_file)
    fprintf(file,"[%s] matrix: %s\n", prog_name, exe_parms->matrix_file);
  fprintf(file,"[%s] matrix: n = (%d)\tm = (%d)\t\tinput_nnz = (%ld)\tnnz_final=(%ld)\n", 
	  prog_name, self->n, self->m, self->nnz_input, self->nnz_final);
  fprintf(file,"[%s] value tolerance (%f)\t\tmarkowitz  tolerance (%f)\tdensity tolerance(%f)\n",
	  prog_name, exe_parms->value_tol, exe_parms->marko_tol, exe_parms->density_tolerance);
  fprintf(file,"[%s] initial extra space (%f)\tinitial inbetween extra space (%f)\n",
	  prog_name, exe_parms->extra_space, exe_parms->extra_space_inbetween);
  fprintf(file,"[%s] #previous pivots (%d)\t\t#min pivots per steps (%d)\n",
	  prog_name, exe_parms->nb_previous_pivots, exe_parms->min_pivot_per_steps);
  fprintf(file,"[%s] #threads (%d)\t\t\t#candidates per bloc (%d)\n",
	  prog_name, exe_parms->nb_threads, exe_parms->nb_candidates_per_block);
  fprintf(file,"[%s] \n", prog_name);

  fprintf(file,"[%s] %d independent sets of pivots were found with a total of %d pivots in %f seconds\n", 
	  prog_name, self->nb_steps, self->sparse_pivots, self->timing_facto_sparse / 1e6);
  if ( ! ( self->reason & TP_reason_dense_too_large ) )
    {
      if (self->dense_pivots) {
	fprintf(file,"[%s] %d pivots were handeled by the dense code in %f seconds\n", 
		prog_name, self->dense_pivots, self->timing_facto_dense / 1e6); 
	fprintf(file,"[%s] The convertion of the sparse Schur to a dense matrix of size %d took %f seconds\n",
		prog_name, self->dense_pivots, self->timing_convert_schur / 1e6);
      }
      
      fprintf(file,"[%s] \n", prog_name);
      fprintf(file,"[%s] The solve was performed in %f seconds\n", prog_name, self->timing_solve/1e6);
      fprintf(file,"[%s] The solve on L took %f seconds ", prog_name, self->timing_solve_L/1e6);
      if (self->dense_pivots)  fprintf(file,", the dense part in %f seconds ", self->timing_solve_dense/1e6);
      fprintf(file,"and the U in %f seconds.\n", self->timing_solve_U/1e6);
    }
  else 
    {
      fprintf(file,"[%s] The factorization was not finished  because the dense part was too large.\n", prog_name);
    }

  fprintf(file,"[%s] The switch to dense code was done ", prog_name);
  if (self->reason & TP_reason_unknown)
    fprintf(file,"for uknown reason.\n");
  else if (self->reason & TP_reason_density)
    fprintf(file,"beacuse the schur became too dense.\n");
  else if (self->reason & TP_reason_no_pivots)
    fprintf(file,"because we did not found engough pivots.\n");
  else if (self->reason & TP_reason_because)
    fprintf(file,"just because.\n");

  if(self->backward_error)
  fprintf(file,"[%s] The backward error is (%e) and the forward error is (%e)\n",
 	  prog_name, self->backward_error, self->forward_error);

  TP_verbose_print_steps(self->stats_first, self->parms);
}

void
create_plot_file(char *plot_file, char *data_file, 
		 char *fig_file, char *title,
		 char *ylabel, int nb_graphs)
{
  FILE *file = fopen(plot_file, "w+");
  int i ;  
 fprintf(file, "set terminal eps\n");
  fprintf(file, "set style data  histogram\n");
  fprintf(file, "set style fill solid border\n");
  fprintf(file, "set style histogram rowstacked\n");
  fprintf(file, "set title \"%s\"\n", title);
  fprintf(file, "set boxwidth 0.6\n");
  fprintf(file, "set output \"%s\"\n", fig_file);
  fprintf(file, "set ylabel \"%s\"\n", ylabel);
  fprintf(file, "set datafile missing \"?\"\n");
  fprintf(file, "plot \'%s\' u 1 t columnheader",data_file);
  for(i = 2; i <= nb_graphs; i++)
    fprintf(file, ",  \'\' u %d t columnheader",i);
  fprintf(file, "\n");
  fclose(file);
}


void
create_time_data(char *filename, TP_verbose verbose)
{
  FILE *file = fopen(filename, "w+");
  TP_verbose_per_step step = verbose->stats_first; 

  fprintf(file, "search\tupdate\tdense\n");
  while(step) 
    {
      fprintf(file, "%f\t%f\t?\n", 
	      step->timing_pivot_search,
	      step->timing_apply_perms);
      step = step->next;
    }
  fprintf(file, "?\t?\t%f\n", verbose->timing_facto_dense);

  fclose(file);
}


void
create_pivots_data(char *filename, TP_verbose verbose)
{
  FILE *file = fopen(filename, "w+");
  TP_verbose_per_step step = verbose->stats_first;

  fprintf(file, "\"#pivot per step\"\t\"dense pivots\"\n");
  while(step)
    {
      fprintf(file, "%d\t?\n", step->nb_pivots);
      step = step->next;
    }
  fprintf(file, "?\t%d\n", verbose->dense_pivots);

  fclose(file);
}


void
TP_verbose_create_dirs_V0(char *dir) 
{
  char *new_dir = malloc((size_t) (strlen(dir)+200) * sizeof(*new_dir));

  *new_dir = '\0';
  strcat(new_dir, dir);
  strcat(new_dir, "/data");
  mkdir(new_dir, 0755);

  *new_dir = '\0';
  strcat(new_dir, dir);
  strcat(new_dir, "/plot");
  mkdir(new_dir, 0755);

  *new_dir = '\0';
  strcat(new_dir, dir);
  strcat(new_dir, "/fig");
  mkdir(new_dir, 0755);

  free(new_dir);
}


void
TP_verbose_draw_graph_V0(TP_verbose verbose)
{
  TP_verbose_parms parms = verbose->parms;
  char *output_dir = parms->output_dir;
  char data_file[2048];
  char plot_file[2048];
  char fig_file[2048];

  /* Handle the time graphs */
  snprintf(plot_file, 2048, "%s/plot/%s_time.gplt", output_dir, parms->outfiles_prefix);
  snprintf(data_file, 2048, "../data/%s_time.dat",  parms->outfiles_prefix);
  snprintf(fig_file,  2048, "../fig/%s_time.eps",   parms->outfiles_prefix);

  create_plot_file(plot_file, data_file, fig_file, "time per pivot set", "time", 3);
  snprintf(data_file, 2048, "%s/data/%s_time.dat", output_dir, parms->outfiles_prefix);
  create_time_data(data_file, verbose);

  /* Handle the #pivots graphs */
  snprintf(plot_file, 2048, "%s/plot/%s_pivots.gplt", output_dir, parms->outfiles_prefix);
  snprintf(data_file, 2048, "../data/%s_pivots.dat",  parms->outfiles_prefix);
  snprintf(fig_file,  2048, "../fig/%s_pivots.eps",   parms->outfiles_prefix);

  create_plot_file(plot_file, data_file, fig_file, "#pivots per pivot set", "#pivots", 2);
  snprintf(data_file, 2048, "%s/data/%s_pivots.dat", output_dir, parms->outfiles_prefix);
  create_pivots_data(data_file, verbose);

  TP_verbose_print_raw(verbose);
}


void 
TP_verbose_destroy_V0(TP_verbose self)
{
  TP_verbose_per_step step = self->stats_first;
  TP_verbose_parms   parms = self->parms;

  while(step)
    {
      TP_verbose_per_step tmp = step->next;
      free(step);
      step = tmp;
    }

  if (!parms->user_out_dir)
    free(parms->output_dir);
  if (parms->user_out_file)
    fclose(parms->out_file);
  free(parms->outfiles_prefix);

  free(parms);
  free(self);
}
