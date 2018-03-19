#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include "TP_auxiliary.h"
#include "TP_paje.h"

#define TP_PAJE_MAXEVENTS  15000
#define TP_PAJE_MAXTPES  12

static char *labels[TP_PAJE_MAXTPES] = {"Get singeltons", "Get eligeble pivtos", "Auxiliary", "Assign score", 
					"First pass", "Discarding Pivots", "Getting Pivots", "Update L",
					"Update S cols", "Update S rows", "Matrix convertion",
					"Dense Factorization"}; 

static double colors[TP_PAJE_MAXTPES][3] = {{1, 1, 0}, {1, 0.8, 0}, {0.9, 0.6, 0}, {0.8, 0.4, 0.15},
					    {0.5, 0.8, 0.2}, {0.2, 0.9, 0.45}, {0, 0, 1}, {0, 1, 0},
					    {0, 1, 1}, {0, 0.8, 0.2}, {0, 0.2, 0.8},  {0.2, 0.4, 1}}; 

typedef struct _TP_paje_event TP_paje_event;

struct _TP_paje_event {
  int id_event; 
  double start;
  double end;
};

struct _TP_paje {
  int nb_threads;
  double time_zero;

  int *pending;
  int *nb_events;
  TP_paje_event **events;
};


TP_paje
TP_paje_create(int nb_threads)
{
  TP_paje self = calloc((size_t) 1, sizeof(*self));
  struct timeval time_s; 
  int i;

  gettimeofday(&time_s, NULL);
  self->time_zero = (double) (time_s.tv_sec * 1e6 + time_s.tv_usec);

  self->nb_threads = nb_threads;

    
  self->pending   = calloc((size_t) nb_threads, sizeof(*self->pending));
  self->nb_events = calloc((size_t) nb_threads, sizeof(*self->nb_events));
  self->events = calloc((size_t) nb_threads, sizeof(*self->events));
  for( i = 0; i < self->nb_threads; i++)
    self->events[i] = calloc((size_t) TP_PAJE_MAXEVENTS, sizeof(**self->events));
  
  return self;
}

void
TP_paje_start_event(TP_paje self, int id)
{
  int thread =  omp_get_thread_num();
  TP_paje_event *event = &self->events[thread][self->nb_events[thread]];
  struct timeval time; 
  gettimeofday(&time, NULL);

  if (self->pending[thread])
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Nested events are not supported");
    
  self->pending[thread] = 1; 
  event->start = (double) (time.tv_sec * 1e6 + time.tv_usec) - self->time_zero;
  event->id_event = id;
}

void
TP_paje_stop_event(TP_paje self)
{
  int thread =  omp_get_thread_num();
  TP_paje_event *event = &self->events[thread][self->nb_events[thread]++];
  struct timeval time; 
  gettimeofday(&time, NULL);

  if (!self->pending[thread])
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Event is asked to stop, but is not started yet");

  self->pending[thread] = 0; 
  event->end = (double) (time.tv_sec * 1e6 + time.tv_usec) - self->time_zero;
}

void
TP_paje_create_file(TP_paje self, char *dir, char *prefix)
{
  FILE *file;
  char filename[2048];
  double end_time = 0;
  int i, j;

  snprintf(filename, 2048, "%s/trace/%s.trace", dir, prefix);
  file = fopen(filename, "w+");

  fprintf(file,"%%EventDef PajeDefineContainerType 1\n");
  fprintf(file,"%% Alias string \n");
  fprintf(file,"%% ContainerType string\n"); 
  fprintf(file,"%% Name string \n");
  fprintf(file,"%%EndEventDef\n");
  fprintf(file,"%%EventDef PajeDefineStateType 3\n");
  fprintf(file,"%% Alias string \n");
  fprintf(file,"%% ContainerType string\n"); 
  fprintf(file,"%% Name string \n");
  fprintf(file,"%%EndEventDef \n");
  fprintf(file,"%%EventDef PajeDefineEntityValue 6\n");
  fprintf(file,"%%Alias string  \n");
  fprintf(file,"%% EntityType string  \n");
  fprintf(file,"%% Name string  \n");
  fprintf(file,"%% Color color \n");
  fprintf(file,"%%EndEventDef  \n");
  fprintf(file,"%%EventDef PajeCreateContainer 7\n");
  fprintf(file,"%% Time date  \n");
  fprintf(file,"%% Alias string  \n");
  fprintf(file,"%% Type string  \n");
  fprintf(file,"%% Container string  \n");
  fprintf(file,"%% Name string  \n");
  fprintf(file,"%%EndEventDef  \n");
  fprintf(file,"%%EventDef PajeDestroyContainer 8\n");
  fprintf(file,"%% Time date  \n");
  fprintf(file,"%% Name string  \n");
  fprintf(file,"%% Type string  \n");
  fprintf(file,"%%EndEventDef  \n");
  fprintf(file,"%%EventDef PajeSetState 10\n");
  fprintf(file,"%% Time date  \n");
  fprintf(file,"%% Type string  \n");
  fprintf(file,"%% Container string  \n");
  fprintf(file,"%% Value string  \n");
  fprintf(file,"%%EndEventDef \n");
  fprintf(file,"1 CT_Prog   0       \"Program\"\n");
  fprintf(file,"1 CT_Thread CT_Prog \"Thread\"\n");
  fprintf(file,"3 ST_ThreadState CT_Thread \"Thread State\"\n");

  for( i = 0; i < TP_PAJE_MAXTPES; i++)
    fprintf(file,"6 \"%s\" ST_ThreadState \"%s\" \"%f %f %f\" \n",
	    labels[i], labels[i], colors[i][0], colors[i][1], colors[i][2]);

  fprintf(file,"6 idle ST_ThreadState \"Idle\"  \"1 0 0\"\n");
  fprintf(file,"7 0.000000 C_Prog CT_Prog 0 \"Program\"\n");

  for( i = 0; i < self->nb_threads; i++)
    fprintf(file,"7  0.000000 C_Thread%d CT_Thread C_Prog \"Thread %d\"\n", i, i);

  for( i = 0; i < self->nb_threads; i++)
    for(j = 0; j < self->nb_events[i]; j++) 
      {
	TP_paje_event *event = &self->events[i][j];
	fprintf(file,"10 %f ST_ThreadState C_Thread%d \"%s\"\n", event->start / 1e6, i, labels[event->id_event]);
	fprintf(file,"10 %f ST_ThreadState C_Thread%d idle\n", event->end / 1e6, i);
      }
  
  for( i = 0; i < self->nb_threads; i++)
    end_time = self->events[i][self->nb_events[i]].end > end_time ? self->events[i][self->nb_events[i]].end : end_time;
  
  for( i = 0; i < self->nb_threads; i++)
    fprintf(file,"8 %f C_Thread%d CT_Thread\n", end_time / 1e6, i);
  
  fprintf(file,"%f C_Prog CT_Prog\n", end_time / 1e6);
 										 
  fclose(file);
}


void
TP_paje_destroy(TP_paje self)
{
  int i;
  free(self->pending);
  free(self->nb_events);

  for( i = 0; i < self->nb_threads; i++)
    free(self->events[i]);
  free(self->events);

  free(self);
}
