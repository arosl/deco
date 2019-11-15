#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_HDR       1024
#define TISSUE_COUNT    16
#define P_H2O            0.0567
#define P_STD            1.01325
#define BOTTOM_GAS_IDX   0
#define SURFACE_STOP    -1
#define ASCENT_RATE     10.0
#define M_PER_BAR       10.0
#define GF_MIN           0.3
#define GF_MAX           1.0


typedef struct {
  int      tissue_count;

  double*  n2_a;

  double*  n2_b;

  double*  n2_hl;

  double*  he_a;

  double*  he_b;

  double*  he_hl;
} par_t;


typedef struct {
  double*  p_n2;

  double*  p_he;

  double   time;

  double*  gas_use;
} state_t;


typedef struct {
  double  f_n2;

  double  f_he;

  double  max_depth;

  int     is_optional;
} gas_t;


typedef struct {
  int      fixed_stop_count;

  double*  fixed_stops;

  double   stop_step;

  int      gas_count;

  gas_t**  gasses;

  double   p_amb;

  double   gf_min;

  double   gf_max;

  int      smooth_gradient;

  double   ascent_rate;

  int      is_derived;

  double   gas_use_dive;

  double   gas_use_deco;

  double   init_time;

  double   init_depth;

  int      init_gas;

  double   finish_time;

  double   finish_depth;

  int      finish_gas;
} schedule_t;


typedef struct {
  void  (*start_table)( int  cols );

  void  (*end_table)( void );

  void  (*start_row)( void );

  void  (*end_row)( void );

  void  (*write_num)( double  d );

  void  (*write_str)( char*  s,
                      int    width );

  void  (*rule)( void );

  void  (*drule)( void );
} printer_t;


typedef struct {
  schedule_t*  schedule;

  double       depth;

  int          time_cnt;

  double*      times;

  printer_t*   printer;

  int          total_is_deco;
} context_t;


typedef struct {
  int      stop_count;

  double*  stop_times;

  double   total_time;

  double*  gas_use;
} deco_t;


typedef struct {
  int      count;

  double*  depths;

  int*     gasses;

  int      finish_idx;
} stop_list_t;


typedef struct {
  int           time_cnt;

  stop_list_t*  stop_list;

  deco_t**      decos;
} table_t;


par_t*  par_create( void )
{
  par_t*  par  =  (par_t*) malloc(sizeof(par_t));
  int     i;

  double n2_a[TISSUE_COUNT]   =  { 1.1696, 1.0000, 0.8618, 0.7562, 0.6667, 0.5600, 0.4947, 0.4500,
                                   0.4187, 0.3798, 0.3497, 0.3223, 0.2850, 0.2737, 0.2523, 0.2327 };
  double n2_b[TISSUE_COUNT]   =  { 0.5578, 0.6514, 0.7222, 0.7825, 0.8126, 0.8434, 0.8693, 0.8910,
                                   0.9092, 0.9222, 0.9319, 0.9403, 0.9477, 0.9544, 0.9602, 0.9653 };
  double n2_hl[TISSUE_COUNT]  =  {    5.0,    8.0,   12.5,   18.5,   27.0,   38.3,   54.3,   77.0,
                                    109.0,  146.0,  187.0,  239.0,  305.0,  390.0,  498.0,  635.0 };
  double he_a[TISSUE_COUNT]   =  { 1.6189, 1.3830, 1.1919, 1.0458, 0.9220, 0.8205, 0.7305, 0.6502,
                                   0.5950, 0.5545, 0.5333, 0.5189, 0.5181, 0.5176, 0.5172, 0.5119 };
  double he_b[TISSUE_COUNT]   =  { 0.4770, 0.5747, 0.6527, 0.7223, 0.7582, 0.7957, 0.8279, 0.8553,
                                   0.8757, 0.8903, 0.8997, 0.9073, 0.9122, 0.9171, 0.9217, 0.9267 };
  double he_hl[TISSUE_COUNT]  =  {   1.88,   3.02,   4.72,   6.99,  10.21,  14.48,  20.53,  29.11,
                                    41.20,  55.19,  70.69,  90.34, 115.29, 147.42, 188.24, 240.03 };

  par->tissue_count  =  TISSUE_COUNT;
  par->n2_a          =  (double*) malloc(TISSUE_COUNT * sizeof(double));
  par->n2_b          =  (double*) malloc(TISSUE_COUNT * sizeof(double));
  par->n2_hl         =  (double*) malloc(TISSUE_COUNT * sizeof(double));
  par->he_a          =  (double*) malloc(TISSUE_COUNT * sizeof(double));
  par->he_b          =  (double*) malloc(TISSUE_COUNT * sizeof(double));
  par->he_hl         =  (double*) malloc(TISSUE_COUNT * sizeof(double));

  for (i = 0; i < TISSUE_COUNT; i++) {
    par->n2_a[i]   =  n2_a[i];
    par->n2_b[i]   =  n2_b[i];
    par->n2_hl[i]  =  n2_hl[i];
    par->he_a[i]   =  he_a[i];
    par->he_b[i]   =  he_b[i];
    par->he_hl[i]  =  he_hl[i];
  }

  return  par;
}


state_t*  state_create( schedule_t*  schedule,
                        par_t*       par )
{
  state_t*  state  =  (state_t*) malloc(sizeof(state_t));

  state->p_n2     =  (double*) malloc(par->tissue_count * sizeof(double));
  state->p_he     =  (double*) malloc(par->tissue_count * sizeof(double));
  state->gas_use  =  (double*) malloc(schedule->gas_count * sizeof(double));

  return state;
}


void  state_free( state_t*  state )
{
  free(state->p_n2);
  free(state->p_he);
  free(state->gas_use);
  free(state);
}


void  state_copy( schedule_t*  schedule,
                  par_t*       par,
                  state_t*     dst,
                  state_t*     src )
{
  int  i;

  for (i = 0; i < par->tissue_count; i++) {
    dst->p_n2[i]  =  src->p_n2[i];
    dst->p_he[i]  =  src->p_he[i];
  }

  for (i = 0; i < schedule->gas_count; i++) {
    dst->gas_use[i]  =  src->gas_use[i];
  }

  dst->time = src->time;
}


void  state_print( schedule_t*  schedule,
                   par_t*       par,
                   state_t*     state )
{
  int  i;

  for (i = 0; i < par->tissue_count; i++) {
    printf(" %5.2f", state->p_n2[i]);
  }
  printf("\n");

  printf(" %5.2f", state->gas_use[0] * schedule->gas_use_dive);
  for (i = 1; i < schedule->gas_count; i++) {
    printf(" %5.2f", state->gas_use[i] * schedule->gas_use_deco);
  }
  printf("\n");
}


void  par_free( par_t*  par )
{
  free(par->n2_a);
  free(par->n2_b);
  free(par->n2_hl);
  free(par->he_a);
  free(par->he_b);
  free(par->he_hl);
  free(par);
}


gas_t*  gas_create( double  f_n2,
                    double  f_he,
                    double  max_depth )
{
  gas_t*  gas  =  (gas_t*) malloc(sizeof(gas_t));

  gas->f_n2       =  f_n2;
  gas->f_he       =  f_he;
  gas->max_depth  =  max_depth;

  return  gas;
}


void  gas_free( gas_t*  gas )
{
  free(gas);
}


double  gas_get_o2( gas_t*  gas )
{
  return  1.0 - gas->f_n2 - gas->f_he;
}


void  trim_gas_str( char*  str )
{
  int  i;

  for (i = 0; str[i] != '\0' && str[i] != ','; i++) {
    ;
  }

  str[i] = '\0';
}


gas_t*  gas_parse( char*  str,
                   int    is_deco,
                   int*   ptr )
{
  gas_t*  gas          =  NULL;
  double  max_o2       =  1.4;
  double  o2;
  double  he;
  double  d;
  int     is_optional  =  0;
  int     n2;
  int     n            =  0;
  int     p            =  0;

  if (ptr != NULL) {
    p = *ptr;
  }

  if (is_deco) {
    max_o2 = 1.6;
  }

  if (str[p] == '[') {
    is_optional = 1;
    (*ptr)++;
    p++;
  }

  if (strncmp(str + p, "air", 3) == 0) {
    n    =  3;
    gas  =  gas_create(0.79, 0.00, 36.0);
  }
  else if (strncmp(str + p, "o2", 2) == 0) {
    n    =  2;
    gas  =  gas_create(0.00, 0.00, 6);
  }
  else if (sscanf(str + p, "ean%lf%n", &o2, &n) == 1) {
    double  max_depth;

    o2         /=  100.0;
    max_depth   =  (max_o2 / o2 - 1.0) * 10.0;

    if (!is_deco && max_depth > 30) {
      max_depth = 30;
    }
    if (is_deco && max_depth > 36) {
      max_depth = 36;
    }

    gas = gas_create(1.0 - o2, 0.00, max_depth);
  }
  else if (sscanf(str + p, "%lf/%lf%n", &o2, &he, &n) == 2) {
    double  max_depth;
    double  n2;

    o2         /=  100.0;
    he         /=  100.0;
    n2          =  1.0 - o2 - he;
    max_depth   =  (max_o2 / o2 - 1.0) * 10.0;

    if (n2 > 0.0) {
      double  n2_depth  =  (3.16 / n2 - 1.0) * 10.0;

      if (n2_depth < max_depth) {
        max_depth = n2_depth;
      }
    }

    gas = gas_create(n2, he, max_depth);
  }

  if (gas != NULL && sscanf(str + p + n, "@%lf%n", &d, &n2) == 1) {
    n               +=  n2;
    gas->max_depth  =   d;
  }

  if (gas == NULL ||
      (!is_optional && str[p + n] != '\0' && str[p + n] != ',') ||
      (is_optional && str[p + n] != ']') ||
      gas->f_n2 < 0.0 || gas->f_he < 0.0 || gas_get_o2(gas) < 0.0) {
    trim_gas_str(str + p);
    fprintf(stderr, "Unknown gas: %s\n", str + p);
    exit(1);
  }

  if (is_optional) {
    n++;
  }

  if (str[p + n] == ',') {
    n++;
  }

  if (ptr != NULL) {
    *ptr += n;
  }

  gas->is_optional  =  is_optional;

  return  gas;
}


int  print_gas_frac( char*   str,
                     double  f )
{
  double  frac;

  f     *=  100.0;
  frac   =  f;
  frac  -=  floor(frac);

  if (frac < 0.05 || frac > 0.95) {
    return  sprintf(str, "%.f", f);
  }
  else {
    return  sprintf(str, "%.1f", f);
  }
}


int  gas_print( char*   str,
                gas_t*  gas )
{
  double  f_o2  =  gas_get_o2(gas);

  if (gas->f_he == 0.0) {
    if (gas->f_n2 >= 0.789 && gas->f_n2 <= 0.791) {
      return sprintf(str, "Air");
    }
    else if (gas_get_o2(gas) >= 0.99) {
      return sprintf(str, "O2");
    }
    else {
      int  n = sprintf(str, "EAN");

      return  print_gas_frac(str + n, f_o2) + n;
    }
  }
  else {
    int  n  =  print_gas_frac(str, f_o2);

    n += sprintf(str + n, "/");

    return print_gas_frac(str + n, gas->f_he) + n;
  }
}


schedule_t*  schedule_create( int     fixed_stop_count,
                              int     gas_count,
                              double  p_amb,
                              double  gf_min,
                              double  gf_max,
                              int     smooth_gradient,
                              double  ascent_rate )
{
  schedule_t*  schedule  =  (schedule_t*) malloc(sizeof(schedule_t));

  schedule->fixed_stop_count  =  fixed_stop_count;
  if (fixed_stop_count == 0) {
    schedule->fixed_stops     =  NULL;
  }
  else {
    schedule->fixed_stops     =  (double*) malloc(fixed_stop_count * sizeof(double));
  }
  schedule->gas_count         =  gas_count;
  schedule->gasses            =  (gas_t**) malloc(gas_count * sizeof(gas_t*));
  schedule->p_amb             =  p_amb;
  schedule->gf_min            =  gf_min;
  schedule->gf_max            =  gf_max;
  schedule->smooth_gradient   =  smooth_gradient;
  schedule->ascent_rate       =  ascent_rate;
  schedule->is_derived        =  0;

  return  schedule;
}


void  schedule_free( schedule_t*  schedule )
{

  if (!schedule->is_derived) {
    int  i;

    for (i = 0; i < schedule->gas_count; i++) {
      gas_free(schedule->gasses[i]);
    }
    free(schedule->fixed_stops);
  }

  free(schedule->gasses);
  free(schedule);
}


double  schedule_get_stop( schedule_t*  schedule,
                           int          index )
{
  int  fixed_count  =  schedule->fixed_stop_count;

  if (index == SURFACE_STOP) {
    return  0.0;
  }

  if (fixed_count == 0) {
    return  schedule->stop_step * (1 + index);
  }

  if (index < fixed_count) {
    return  schedule->fixed_stops[index];
  }

  return  schedule->fixed_stops[fixed_count - 1] + schedule->stop_step * (index - fixed_count + 1);
}


gas_t*  schedule_get_gas( schedule_t*  schedule,
                          int          index )
{
  if (index == schedule->gas_count) {
    return  NULL;
  }

  return  schedule->gasses[index];
}


int  schedule_get_deepest_stop( schedule_t*  schedule,
                                double       depth )
{
  int  i;

  for (i = 0; schedule_get_stop(schedule, i) < depth; i++) {
    ;
  }

  return  i - 1;
}


schedule_t*  schedule_derive( schedule_t*  schedule,
                              int          index )
{
  schedule_t*  new_schedule;
  int          i;

  if (!schedule->gasses[index]->is_optional) {
    return  NULL;
  }

  new_schedule  =  schedule_create(0,
                                   schedule->gas_count - 1,
                                   schedule->p_amb,
                                   schedule->gf_min,
                                   schedule->gf_max,
                                   schedule->smooth_gradient,
                                   schedule->ascent_rate);

  for (i = 0; i < new_schedule->gas_count; i++) {
    if (i < index) {
      new_schedule->gasses[i] = schedule->gasses[i];
    }
    else {
      new_schedule->gasses[i] = schedule->gasses[i + 1];
    }
  }

  new_schedule->fixed_stop_count  =  schedule->fixed_stop_count;
  new_schedule->fixed_stops       =  schedule->fixed_stops;
  new_schedule->stop_step         =  schedule->stop_step;
  new_schedule->is_derived        =  1;
  new_schedule->gas_use_dive      =  schedule->gas_use_dive;
  new_schedule->gas_use_deco      =  schedule->gas_use_deco;
  new_schedule->init_time         =  schedule->init_time;
  new_schedule->init_depth        =  schedule->init_depth;
  new_schedule->init_gas          =  schedule->init_gas;
  new_schedule->finish_time       =  schedule->finish_time;
  new_schedule->finish_depth      =  schedule->finish_depth;
  new_schedule->finish_gas        =  schedule->finish_gas;

  return  new_schedule;
}


stop_list_t*  stop_list_create( schedule_t*  schedule,
                                double       depth )
{
  stop_list_t*  list      =  (stop_list_t*) malloc(sizeof(stop_list_t));
  int           idx       =  schedule_get_deepest_stop(schedule, depth);
  int           gas_idx   =  BOTTOM_GAS_IDX;
  gas_t*        next_gas  =  schedule_get_gas(schedule, gas_idx + 1);
  int           i;

  list->count       =  idx + 1;
  list->depths      =  (double*) malloc(list->count * sizeof(double));
  list->gasses      =  (int*) malloc(list->count * sizeof(int));
  list->finish_idx  =  -1;

  for (i = list->count - 1; i >= 0; i--) {
    double  depth  =  schedule_get_stop(schedule, i);

    while (next_gas != NULL && depth <= next_gas->max_depth) {
      gas_idx++;
      next_gas  =  schedule_get_gas(schedule, gas_idx + 1);
    }

    list->depths[i] = depth;
    list->gasses[i] = gas_idx;
  }

  if (schedule->finish_time != 0.0) {
    int  do_insert  =  0;

    for (i = list->count - 1; i >= 0; i--) {
      if (list->depths[i] == schedule->finish_depth) {
        list->finish_idx = i;

        if (schedule->finish_gas != -1 && list->gasses[i] != schedule->finish_gas) {
          list->finish_idx++;
          do_insert = 1;
        }

        break;
      }

      if (schedule->finish_depth > list->depths[i]) {
        list->finish_idx = i + 1;

        do_insert = 1;

        break;
      }
    }

    if (list->finish_idx == -1) {
      // finish profile shallower than last stop
      list->finish_idx = 0;

      do_insert = 1;
    }

    if (do_insert) {
      int  j;

      list->count++;
      list->depths  =  (double*) realloc(list->depths, (list->count + 1) * sizeof(double));
      list->gasses  =  (int*) realloc(list->gasses, (list->count + 1) * sizeof(int));

      for (j = list->count - 1; j > list->finish_idx; j--) {
        list->depths[j] = list->depths[j - 1];
        list->gasses[j] = list->gasses[j - 1];
      }

      list->depths[list->finish_idx] = schedule->finish_depth;

      if (schedule->finish_gas == -1) {
        // check if we can use a shallower gas
        int  gas_idx  =  list->gasses[list->finish_idx];

        while (gas_idx < schedule->gas_count - 1) {
          gas_t*  gas  =  schedule_get_gas(schedule, gas_idx + 1);

          if (gas->max_depth < schedule->finish_depth) {
            break;
          }

          gas_idx++;
          printf("switch gas\n");
        }

        list->gasses[list->finish_idx]  =  gas_idx;
      }
      else {
        list->gasses[list->finish_idx] = schedule->finish_gas;
      }
    }
  }

  return  list;
}


void  stop_list_free( stop_list_t*  list )
{
  free(list->gasses);
  free(list->depths);
  free(list);
}


double  stop_list_get_depth( stop_list_t*  list,
                             int           index )
{
  if (index == SURFACE_STOP) {
    return  0.0;
  }

  return  list->depths[index];
}


double  frac_to_pp( double  f,
                    double  p )
{
  return  f * (p - P_H2O);
}


double  depth_to_p( schedule_t*  schedule,
                    double       d )
{
  return  schedule->p_amb + d / M_PER_BAR;
}


void  saturate( schedule_t*  schedule,
                par_t*       par,
                state_t*     state,
                double       p,
                gas_t*       gas )
{
  double  p_n2  =  frac_to_pp(gas->f_n2, p);
  double  p_he  =  frac_to_pp(gas->f_he, p);
  int     i;

  for (i = 0; i < par->tissue_count; i++) {
    state->p_n2[i]  =  p_n2;
    state->p_he[i]  =  p_he;
  }

  for (i = 0; i < schedule->gas_count; i++) {
    state->gas_use[i]  =  0;
  }
  state->time = 0.0;
}


void  dive( schedule_t*  schedule,
            par_t*       par,
            state_t*     state,
            double       p1,
            double       p2,
            double       t,
            int          gas_idx )
{
  if (t == 0.0) {
    return;
  }

  gas_t*  gas    =  schedule->gasses[gas_idx];
  double  p1_n2  =  frac_to_pp(gas->f_n2, p1);
  double  p1_he  =  frac_to_pp(gas->f_he, p1);
  double  p2_n2  =  frac_to_pp(gas->f_n2, p2);
  double  p2_he  =  frac_to_pp(gas->f_he, p2);
  double  n2_r   =  (p2_n2 - p1_n2) / t;
  double  he_r   =  (p2_he - p1_he) / t;
  int     i;

  if (t < 0.0) {
    fprintf(stderr, "Time should be non-negative\n");
    exit(1);
  }


  for (i = 0; i < par->tissue_count; i++) {
    double  n2_k  =  M_LN2 / par->n2_hl[i];
    double  he_k  =  M_LN2 / par->he_hl[i];

    state->p_n2[i]  =  p1_n2 + n2_r * (t - 1 / n2_k) - (p1_n2 - state->p_n2[i] - n2_r / n2_k) * exp(-n2_k * t);
    state->p_he[i]  =  p1_he + he_r * (t - 1 / he_k) - (p1_he - state->p_he[i] - he_r / he_k) * exp(-he_k * t);
  }

  state->time += t;

  state->gas_use[gas_idx] += (p1 + p2) / 2.0 * t;
}


double  get_p_min( par_t*    par,
                   state_t*  state,
                   double    gf )
{
  double  p_min  =  0.0;
  int     i;

  for (i = 0; i < par->tissue_count; i++) {
    double  a  =  (par->n2_a[i] * state->p_n2[i] + par->he_a[i] * state->p_he[i]) / (state->p_n2[i] + state->p_he[i]);
    double  b  =  (par->n2_b[i] * state->p_n2[i] + par->he_b[i] * state->p_he[i]) / (state->p_n2[i] + state->p_he[i]);
    double  p  =  (state->p_n2[i] + state->p_he[i] - a * gf) / (gf / b + 1.0 - gf);

    if (i == 0 || p > p_min) {
      p_min = p;
    }
  }

  return  p_min;
}


deco_t*  deco_create( int  stop_count,
                      int  gas_count )
{
  deco_t*  deco  =  (deco_t*) malloc(sizeof(deco_t));

  deco->stop_count  =  stop_count;
  deco->stop_times  =  (double*) malloc(stop_count * sizeof(double));
  deco->gas_use     =  (double*) malloc(gas_count * sizeof(double));

  return  deco;
}


void  deco_free( deco_t*  deco )
{
  free(deco->stop_times);
  free(deco->gas_use);
  free(deco);
}


double  get_gf( schedule_t*  schedule,
                double       depth,
                double       first_stop )
{
  if (isnan(first_stop)) {
    return  schedule->gf_min;
  }

  if (depth > first_stop) {
    printf("dd %f %f\n", depth, first_stop);
  }

  return  (schedule->gf_max * (first_stop - depth) + schedule->gf_min * depth) / first_stop;
}


deco_t*  get_deco( par_t*        par,
                   state_t*      state,
                   schedule_t*   schedule,
                   double        depth,
                   stop_list_t*  stop_list,
                   state_t*      tmp_state )
{
  double   ascent_rate  =  schedule->ascent_rate;
  deco_t*  deco         =  NULL;
  double   d1           =  depth;
  double   p1           =  depth_to_p(schedule, d1);
  double   d2;
  double   p2;
  double   first_stop   =  NAN;
  int      stop_idx;

  stop_idx  =  stop_list->count - 1;

  if (stop_idx == SURFACE_STOP) {
    fprintf(stderr, "First stop deeper than bottom depth\n");
    exit(1);
  }

  d2   =  stop_list_get_depth(stop_list, stop_idx);
  p2   =  depth_to_p(schedule, d2);

  dive(schedule, par, state, p1, p2, (d1 - d2) / ascent_rate, 0);

  do {
    double  atime;
    int     dt  =  0;

    d1  =  d2;
    p1  =  p2;
    d2  =  stop_list_get_depth(stop_list, stop_idx - 1);
    p2  =  depth_to_p(schedule, d2);

    atime  =  (d1 - d2) / ascent_rate;

    if (isnan(first_stop)) {
      state_copy(schedule, par, tmp_state, state);

      if (stop_idx == stop_list->finish_idx) {
        dive(schedule, par, tmp_state, p1, p1, ceil(schedule->finish_time), stop_list->gasses[stop_idx]);
        dt = ceil(schedule->finish_time);
      }
      dive(schedule, par, tmp_state, p1, p2, atime, stop_list->gasses[stop_idx]);

      if (get_p_min(par, tmp_state, get_gf(schedule, d2, first_stop)) > p2) {
        if (!schedule->smooth_gradient) {
          first_stop = d1;
        }
        else {
          int  steps  =  100;
          int  i;

          for (i = 1; i < steps; i++) {
            double d  =  d1 + (d2 - d1) * i / steps;
            double p  =  depth_to_p(schedule, d);

            state_copy(schedule, par, tmp_state, state);
            dive(schedule, par, tmp_state, p1, p, (d1 - d) / ascent_rate, stop_list->gasses[stop_idx]);

            if (get_p_min(par, tmp_state, get_gf(schedule, d, first_stop)) > p) {
              break;
            }
          }

          first_stop = d1 + (d2 - d1) * i / steps;
        }
      }
      else {
        state_copy(schedule, par, state, tmp_state);
      }
    }

    if (!isnan(first_stop)) {
      int  start_i  =  1;

      if (stop_idx == stop_list->finish_idx) {
        start_i = ceil(schedule->finish_time); // TODO: this will actually include ascent to next step
      }
      else if (schedule->finish_time != 0.0 && schedule->finish_depth == d1) {
        // we do not have to make a stop if we already made a finish profile at this depth
        start_i = 0;
      }

      for (dt = start_i; dt < 5000; dt++) {
        state_copy(schedule, par, tmp_state, state);

        if (dt - atime > 0.0) {
          dive(schedule, par, tmp_state, p1, p1, dt - atime, stop_list->gasses[stop_idx]);
        }

        dive(schedule, par, tmp_state, p1, p2, atime, stop_list->gasses[stop_idx]);

        if (get_p_min(par, tmp_state, get_gf(schedule, d2, first_stop)) <= p2) {
          state_copy(schedule, par, state, tmp_state);
          break;
        }
      }

      if (dt == 5000) {
        fprintf(stderr, "The profile and/or gas selection seems to make deco impossible\n");
        exit(1);
      }
    }

    if (deco == NULL) {
      if (dt > 0) {
        deco = deco_create(stop_idx + 1, schedule->gas_count);
      }
    }

    if (deco != NULL) {
      deco->stop_times[stop_idx] = dt;
    }

    stop_idx--;
  } while (stop_idx != SURFACE_STOP);

  if (deco == NULL) {
    deco = deco_create(0, schedule->gas_count);
  }

  return  deco;
}


void  do_nothing( void )
{
}


void  test_rule( void )
{
  printf("------\n");
}


void  start_table_human( int  cols )
{
}


void  end_table_human( void )
{
  printf("\n");
}


void  end_row_human( void )
{
  printf("\n");
}


void  write_num_human( double  d )
{
  printf(" %4.f", d);
}


void  write_str_human( char*  s,
                       int    width )
{
  int  len;
  int  i;

  if (width == 1) {
    printf(" %4.4s", s);
    return;
  }

  width  =  5 * width - 1;
  len    =  strlen(s);

  printf(" ");

  for (i = 0; i < (width - len + 1) / 2; i++) {
    printf(" ");
  }

  printf("%s", s);

  for (i = 0; i < (width - len) / 2; i++) {
    printf(" ");
  }
}


void  drule_human( void )
{
  printf("--\n");
}


printer_t*  printer_create_human( void )
{
  printer_t*  printer  =  (printer_t*) malloc(sizeof(printer_t));

  printer->start_table  =  start_table_human;
  printer->end_table    =  end_table_human;
  printer->start_row    =  do_nothing;
  printer->end_row      =  end_row_human;
  printer->write_num    =  write_num_human;
  printer->write_str    =  write_str_human;
  printer->rule         =  do_nothing;
  printer->drule        =  drule_human;

  return  printer;
}


void  start_table_ps( int  cols )
{
  printf("%d tablestart\n", cols);
}


void  end_table_ps( void )
{
  printf("tableend\n\n");
}


void  start_row_ps( void )
{
  printf("  rowstart\n");
  printf("   ");
}


void  end_row_ps( void )
{
  printf("\n  rowend\n");
}


void  write_num_ps( double  d )
{
  printf(" (%.f) 1 field", d);
}


void  write_str_ps( char*  s,
                    int    width )
{
  printf(" (%s) %d field", s, width);
}


void  rule_ps( void )
{
  printf("  rule\n");
}


void  drule_ps( void )
{
  printf("drule\n");
}


printer_t*  printer_create_ps( void )
{
  printer_t*  printer  =  (printer_t*) malloc(sizeof(printer_t));

  printer->start_table  =  start_table_ps;
  printer->end_table    =  end_table_ps;
  printer->start_row    =  start_row_ps;
  printer->end_row      =  end_row_ps;
  printer->write_num    =  write_num_ps;
  printer->write_str    =  write_str_ps;
  printer->rule         =  rule_ps;
  printer->drule        =  drule_ps;

  return  printer;
}


void  printer_free( printer_t*  printer )
{
  free(printer);
}


void  table_print( printer_t*    printer,
                   schedule_t*   schedule,
                   stop_list_t*  stop_list,
                   stop_list_t*  orig_stop_list,
                   double        depth,
                   double*       times,
                   deco_t**      decos,
                   int           time_cnt,
                   int           total_is_deco )
{
  int  i;
  int  stop_cnt       =  0;
  int  first_gas      =  0;

  for (i = 0; i < time_cnt; i++) {
    if (decos[i]->stop_count > stop_cnt) {
      stop_cnt = decos[i]->stop_count;
    }
  }

  if (orig_stop_list == NULL) {
    char  hdr[MAX_HDR];
    int   n;

    printer->start_row();

    printer->write_str("Dpt", 1);
    printer->write_str("O2", 1);

    n = sprintf(hdr, "%.fm: ", depth);

    for (i = 0; i < schedule->gas_count; i++) {
      if (i > 0) {
        n += sprintf(hdr + n, ", ");
      }
      n += gas_print(hdr + n, schedule->gasses[i]);
    }
    sprintf(hdr + n, ", GF: %.2f / %.2f", schedule->gf_min, schedule->gf_max);
    printer->write_str(hdr, time_cnt);

    printer->end_row();

    if (schedule->init_time != 0.0) {
      printer->start_row();
      printer->write_num(schedule->init_depth);
      printer->write_num(100.0 * gas_get_o2(schedule->gasses[schedule->init_gas]));

      for (i = 0; i < time_cnt; i++) {
        printer->write_num(schedule->init_time);
      }
      printer->end_row();

      if (schedule->init_gas != 0) {
        printer->rule();
      }
    }

    printer->start_row();
    printer->write_num(depth);
    printer->write_num(100.0 * gas_get_o2(schedule->gasses[0]));

    for (i = 0; i < time_cnt; i++) {
      printer->write_num(times[i]);
    }
    printer->end_row();
  }

  int  last_gas  =  -1;

  if (schedule->init_time != 0.0) {
    last_gas  =  0;
  }

  for (i = stop_cnt - 1; i >= 0; i--) {
    double  depth  =  stop_list->depths[i];
    int     gas    =  stop_list->gasses[i];
    int     j;

    if (orig_stop_list != NULL && gas != orig_stop_list->gasses[i]) {
      // indicate that we should print the rest
      first_gas = gas;
      orig_stop_list = NULL;
    }

    if (orig_stop_list == NULL) {
      for (j = 0; j < time_cnt; j++) {
        deco_t*  deco  =  decos[j];

        if (deco->stop_times[i] != 0) {
          break;
        }
      }

      if (j == time_cnt) {
        // this stop is not used at all
        continue;
      }

      if (gas != last_gas) {
        printer->rule();
      }

      printer->start_row();
      printer->write_num(depth);
      printer->write_num(100.0 * gas_get_o2(schedule->gasses[gas]));

      for (j = 0; j < time_cnt; j++) {
        deco_t*  deco  =  decos[j];

        if (i >= deco->stop_count || deco->stop_times[i] == 0) {
          printer->write_str("", 1);
        }
        else {
          printer->write_num(deco->stop_times[i]);
        }
      }
      printer->end_row();
    }

    last_gas = gas;
  }
  printer->rule();

  printer->start_row();

  if (total_is_deco) {
    printer->write_str("Deco ", 2);
  }
  else {
    printer->write_str("Total ", 2);
  }

  for (i = 0; i < time_cnt; i++) {
    double  total  =  decos[i]->total_time;

    if (total_is_deco) {
      total -= schedule->init_time;
      total -= times[i];
    }

    printer->write_num(total);
  }
  printer->end_row();

  if (schedule->gas_use_dive != 0.0) {
    printer->rule();

    for (i = first_gas; i < schedule->gas_count; i++) {
      double  gas_use;
      int     j;

      if (i == 0) {
        gas_use  =  schedule->gas_use_dive;
      }
      else {
        gas_use  =  schedule->gas_use_deco;
      }

      printer->start_row();
      if (i == first_gas) {
        printer->write_str("Use", 1);
      }
      else {
        printer->write_str("", 1);
      }
      printer->write_num(100.0 * gas_get_o2(schedule->gasses[i]));

      for (j = 0; j < time_cnt; j++) {
        printer->write_num(decos[j]->gas_use[i] * gas_use);
      }

      printer->end_row();
    }
  }

}


void  usage( void )
{
  printf("\n");
  printf("usage: deco <options>\n");
  printf("\n");
  printf("Makes a deco table for a specific depth with several bottom times. The model is\n");
  printf("ZHL-16B with gradient factors. The descent is instantaneous and the ascent to\n");
  printf("the first stop is not included in the stop time, but is included in the total\n");
  printf("deco time. Each stop time includes the ascent to the next stop (or the\n");
  printf("surface).\n");
  printf("\n");
  printf("The options are:\n");
  printf("\n");
  printf("  -h Print this message.\n");
  printf("\n");
  printf("  -d Depth in meters. (Required)\n");
  printf("\n");
  printf("  -t Comma separated bottom times. (Required)\n");
  printf("\n");
  printf("  -g Comma separated gasses, first bottom gas, then deco gasses. Some examples:\n");
  printf("     air,o2,ean32,18/45. A deco gas surrounded by square brackets is optional\n");
  printf("     resulting in a table variant without that gas. (Default air)\n");
  printf("\n");
  printf("  -f Comma separated gradient factors. Append an 'r' for ragged start of\n");
  printf("     gradient which means that the low gradient factor is always on a deco\n");
  printf("     stop. The default is that the gradient can start anywhere. (Default\n");
  printf("     %.1f,%.1f)\n", GF_MIN, GF_MAX);
  printf("\n");
  printf("  -s Comma separated stop depths in meters starting with the shallowest. If the\n");
  printf("     last depth is less than or equal to the second to last, it is taken to be\n");
  printf("     the step size for the rest of the stops. If not, the step size is 3\n");
  printf("     meters. Some examples: 3 (3, 6, 9, ...), 5,10,2 (5, 10, 12, 14, ...), 6\n");
  printf("     (6, 9, 12, ...), 1,1 (1, 2, 3, ...). (Default 3)\n");
  printf("\n");
  printf("  -a Ambient pressure in bars. (Default %.5f)\n", P_STD);
  printf("\n");
  printf("  -r Ascent rate in m/min. (Default %.1f)\n", ASCENT_RATE);
  printf("\n");
  printf("  -u Comma separated gas usage for the dive and the deco. With this option, gas\n");
  printf("     usage is calculated. If only one number is given, it is applied to both\n");
  printf("     dive and deco.\n");
  printf("\n");
  printf("  -i Initial profile give as depth comma time. This is a square profile\n");
  printf("     inserted at the start of the dive followed by an instantaneous move to the\n");
  printf("     given dive depth and the rest of the dive.\n");
  printf("\n");
  printf("  -j Finish profile given as depth comma time. This works like a deco stop of\n");
  printf("     the given minimum time. If deeper than the bottom depth, the depth just\n");
  printf("     switches to that depth after the bottom time.\n");
  printf("\n");
  printf("  -k Index of gas used in initial profile. The index refers to the gas list and\n");
  printf("     starts from zero. (Default 0, bottom gas)\n");
  printf("\n");
  printf("  -l Index of gas used in finish profile. The index refers to the gas list and\n");
  printf("     starts from zero. (Default is current deco gas if applicable, otherwise\n");
  printf("     bottom gas)\n");
  printf("\n");
  printf("  -z Count the whole dive instead of just the deco as total time.\n");
  printf("\n");
  printf("  -p Output in PostScript format (appropriate header and footer needed).\n");
  printf("\n");
  printf("  --hdr Output in PostScript header (no other options allowed).\n");
  printf("\n");
  printf("  --jnr Output in PostScript joiner for putting two tables on the same page (no\n");
  printf("        other options allowed).\n");
  printf("\n");
  printf("  --ftr Output in PostScript footer (no other options allowed).\n");
  printf("\n");
  printf("Example:\n");
  printf("\n");
  printf("  deco -d 60 -t 15,20,25,30,40,50,60 -g 18/55,ean32,o2 -f 0.3,0.9 -s 6\n");
  printf("\n");
  printf("This creates a table for dives to 60 meters using trimix 18/55 for the bottom\n");
  printf("gas and EAN32 and O2 for deco. The last deco stop is 6 meters.\n");
  printf("\n");

  exit(1);
}


void  print_header( void )
{
  printf("%%!PS-Adobe-2.0\n");
  printf("%%BoundingBox: 0 0 595 842\n");
  printf("\n");
  printf("595 0 translate\n");
  printf("90 rotate\n");
  printf("\n");
  printf("0.7 setlinewidth\n");
  printf("\n");
  printf("/Helvetica-bold findfont\n");
  printf("6.5 scalefont\n");
  printf("setfont\n");
  printf("\n");
  printf("/halfheight\n");
  printf("  newpath 0 0 moveto (M) true charpath flattenpath pathbbox\n");
  printf("  exch pop exch sub exch pop 2 div newpath\n");
  printf("def\n");
  printf("\n");
  printf("/cshow {\n");
  printf("  dup stringwidth pop 2 div neg halfheight neg rmoveto show\n");
  printf("} def\n");
  printf("\n");
  printf("/lshow {\n");
  printf("  dup stringwidth pop neg halfheight neg rmoveto show\n");
  printf("} def\n");
  printf("\n");
  printf("/rshow {\n");
  printf("  0 halfheight neg rmoveto show\n");
  printf("} def\n");
  printf("\n");
  printf("/xbase 150 def\n");
  printf("\n");
  printf("/ybase 500 def\n");
  printf("\n");
  printf("/xdiff 223 def\n");
  printf("/ydiff 390 def\n");
  printf("\n");
  printf("/yst ydiff 60 div def\n");
  printf("\n");
  printf("/ygap 1.5 def\n");
  printf("\n");
  printf("/x xbase def\n");
  printf("\n");
  printf("/y ybase def\n");
  printf("\n");
  printf("/heading {\n");
  printf("    xbase xdiff 2 div add ybase 20 add moveto\n");
  printf("    gsave\n");
  printf("    1.6 dup scale\n");
  printf("    cshow\n");
  printf("    grestore\n");
  printf("} def\n");
  printf("\n");
  printf("/drawlines {\n");
  printf("    /y1 sy yst 2 div add def\n");
  printf("    /y2 y yst 2 div add def\n");
  printf("\n");
  printf("    xbase y1 moveto\n");
  printf("    xbase y2 lineto\n");
  printf("\n");
  printf("    xbase xdiff add y2 lineto\n");
  printf("    xbase xdiff add y1 lineto\n");
  printf("\n");
  printf("    closepath\n");
  printf("\n");
  printf("    xbase xst add y1 moveto\n");
  printf("    xbase xst add tot yst 2 div add lineto\n");
  printf("\n");
  printf("    xbase xst add y2 moveto\n");
  printf("    xbase xst add tot yst 2 div sub lineto\n");
  printf("\n");
  printf("    xbase xst 2 mul add y1 moveto\n");
  printf("    xbase xst 2 mul add y2 lineto\n");
  printf("\n");
  printf("    %% add 1 for rounding issues\n");
  printf("    xst 5 mul xst 3 mul xst colcnt 1 sub mul 1 add {\n");
  printf("        xbase add dup\n");
  printf("        y1 isfull {yst sub} if moveto\n");
  printf("        y2 lineto\n");
  printf("    } for\n");
  printf("\n");
  printf("    isfull {\n");
  printf("        xbase y1 yst sub moveto\n");
  printf("        xbase xst colcnt mul add y1 yst sub lineto\n");
  printf("    } if\n");
  printf("\n");
  printf("    stroke\n");
  printf("} def\n");
  printf("\n");
  printf("/tablestart {\n");
  printf("    dup /colcnt exch def\n");
  printf("    xdiff exch div /xst exch def\n");
  printf("\n");
  printf("    /isfull true def\n");
  printf("    /sy y def\n");
  printf("} def\n");
  printf("\n");
  printf("/tableend {\n");
  printf("    drawlines\n");
  printf("    /y y yst sub ygap sub def\n");
  printf("} def\n");
  printf("\n");
  printf("/drule {\n");
  printf("    drawlines\n");
  printf("    /isfull false def\n");
  printf("    /y y ygap sub def\n");
  printf("    /sy y def\n");
  printf("} def\n");
  printf("\n");
  printf("/rowstart {\n");
  printf("\n");
  printf("} def\n");
  printf("\n");
  printf("/rowend {\n");
  printf("    /y y yst sub def\n");
  printf("    /x xbase def\n");
  printf("} def\n");
  printf("\n");
  printf("/field {\n");
  printf("    dup 2 eq {/tot y def} if\n");
  printf("    x y moveto\n");
  printf("    dup xst mul x add /x exch def\n");
  printf("\n");
  printf("    xst mul 2 div 0 rmoveto\n");
  printf("    cshow\n");
  printf("} def\n");
  printf("\n");
  printf("/rule {\n");
  printf("    xbase y yst 2 div add moveto\n");
  printf("    xbase xst colcnt mul add y yst 2 div add lineto\n");
  printf("} def\n");
  printf("\n");
  printf("\n");
  printf("\n");
  printf("/margin 5 def\n");
  printf("\n");
  printf("/corner {\n");
  printf("    margin neg margin rmoveto\n");
  printf("\n");
  printf("    save\n");
  printf("    -5 0 rmoveto\n");
  printf("    -10 0 rlineto\n");
  printf("    stroke\n");
  printf("    restore\n");
  printf("\n");
  printf("    0 5 rmoveto\n");
  printf("    0 10 rlineto\n");
  printf("    stroke\n");
  printf("} def\n");
  printf("\n");
  printf("/halfcorner {\n");
  printf("    margin neg margin rmoveto\n");
  printf("\n");
  printf("    0 5 rmoveto\n");
  printf("    0 10 rlineto\n");
  printf("    stroke\n");
  printf("} def\n");
  printf("\n");
  printf("/halfcorner2 {\n");
  printf("    margin neg margin rmoveto\n");
  printf("\n");
  printf("    -5 0 rmoveto\n");
  printf("    -10 0 rlineto\n");
  printf("    stroke\n");
  printf("} def\n");
  printf("\n");
  printf("/corners {\n");
  printf("    %% add 0.1 to avoid rounding errors\n");
  printf("    ybase ydiff sub y yst add ygap add 0.1 add gt {\n");
  printf("        xbase xdiff 2 div add ybase yst add moveto\n");
  printf("        (OVERFLOW) cshow\n");
  printf("    } if\n");
  printf("\n");
  printf("    save\n");
  printf("    0.35 setlinewidth\n");
  printf("\n");
  printf("    xbase ybase 30 add yst 2 div add moveto\n");
  printf("    corner\n");
  printf("\n");
  printf("    xbase 40 add ybase 30 add yst 2 div add moveto\n");
  printf("    halfcorner\n");
  printf("\n");
  printf("    xbase ybase yst 2 div add moveto\n");
  printf("    halfcorner2\n");
  printf("\n");
  printf("    xbase xdiff add ybase yst 2 div add moveto\n");
  printf("    save\n");
  printf("    -90 rotate\n");
  printf("    halfcorner\n");
  printf("    restore\n");
  printf("\n");
  printf("    xbase xdiff add ybase 30 add yst 2 div add moveto\n");
  printf("    save\n");
  printf("    -90 rotate\n");
  printf("    corner\n");
  printf("    restore\n");
  printf("\n");
  printf("    xbase xdiff add 40 sub ybase 30 add yst 2 div add moveto\n");
  printf("    save\n");
  printf("    -90 rotate\n");
  printf("    halfcorner2\n");
  printf("    restore\n");
  printf("\n");
  printf("    xbase ybase yst 2 div add ydiff sub moveto\n");
  printf("    save\n");
  printf("    90 rotate\n");
  printf("    corner\n");
  printf("    restore\n");
  printf("\n");
  printf("    xbase 40 add ybase yst 2 div add ydiff sub moveto\n");
  printf("    save\n");
  printf("    90 rotate\n");
  printf("    halfcorner2\n");
  printf("    restore\n");
  printf("\n");
  printf("    xbase xdiff add 40 sub ybase yst 2 div add ydiff sub moveto\n");
  printf("    save\n");
  printf("    180 rotate\n");
  printf("    halfcorner\n");
  printf("    restore\n");
  printf("\n");
  printf("    xbase xdiff add ybase yst 2 div add ydiff sub moveto\n");
  printf("    save\n");
  printf("    180 rotate\n");
  printf("    corner\n");
  printf("    restore\n");
  printf("    restore\n");
  printf("} def\n");
}


void  print_joiner( void )
{
  printf("\n");
  printf("corners\n");
  printf("\n");
  printf("xdiff margin 2 mul add 25 add 0 translate\n");
  printf("\n");
  printf("/x xbase def\n");
  printf("\n");
  printf("/y ybase def\n");
  printf("\n");
}


void  print_footer( void )
{
  printf("\n");
  printf("corners\n");
  printf("\n");
  printf("showpage\n");
}


double*  parse_double_list( char*  s,
                            int*   count_ptr )
{
  double*  array;
  int      p    =  0;
  int      cnt  =  0;
  int      n;

  while (s[p] != '\0') {
    double  dummy;

    if (sscanf(s + p, p == 0 ? "%lf%n" : ",%lf%n", &dummy, &n) != 1 || dummy <= 0.0) {
      return  NULL;
    }

    cnt++;
    p += n;
  }

  array  =  (double*) malloc(cnt * sizeof(double));
  cnt    =  0;
  p      =  0;

  while (s[p] != '\0') {
    sscanf(s + p, p == 0 ? "%lf%n" : ",%lf%n", array + cnt, &n);
    cnt++;
    p += n;
  }

  *count_ptr  =  cnt;

  return  array;
}


gas_t**  parse_gas_list( char*  s,
                         int*   count_ptr )
{
  gas_t**  array;
  int      p        =  0;
  int      cnt      =  0;
  int      is_deco  =  0;

  while (s[p] != '\0') {
    free(gas_parse(s, 1, &p));
    cnt++;
  }

  array  =  (gas_t**) malloc(cnt * sizeof(gas_t*));
  cnt    =  0;
  p      =  0;

  while (s[p] != '\0') {
    array[cnt] = gas_parse(s, is_deco, &p);
    is_deco = 1;
    cnt++;
  }

  *count_ptr  =  cnt;

  return  array;
}


context_t*  parse_options( int      argc,
                           char**   argv )
{
  schedule_t*  schedule;
  context_t*   context;
  double       depth            =  0.0;
  int          gas_cnt          =  0;
  gas_t**      gasses           =  NULL;
  double       gf_min           =  0.0;
  double       gf_max           =  0.0;
  double       smooth_gradient  =  1;
  double       p_amb            =  0.0;
  double       arate            =  0.0;
  int          time_cnt         =  0;
  double*      times            =  NULL;
  int          input_stop_cnt   =  0;
  double*      input_stops      =  NULL;
  double       stop_step        =  3.0;
  double       gas_use_dive     =  0.0;
  double       gas_use_deco     =  0.0;
  printer_t*   printer          =  NULL;
  double       init_depth       =  0.0;
  double       init_time        =  0.0;
  int          init_gas         =  -1;
  double       finish_depth     =  0.0;
  double       finish_time      =  0.0;
  int          finish_gas       =  -1;
  int          total_is_deco    =  1;
  int          i;
  int          n;

  if (argc == 1) {
    usage();
  }

  for (i = 1; i < argc; i++) {
    char*  s  =  argv[i];
    char*  o  =  argv[i + 1]; // only used when i + 1 < argc

    if (strcmp(s, "-h") == 0) {
      usage();
    }
    else if (strncmp(s, "--", 2) == 0) {
      if (argc != 2) {
        usage();
      }

      if (strcmp(s, "--hdr") == 0) {
        print_header();
      }
      else if (strcmp(s, "--jnr") == 0) {
        print_joiner();
      }
      else if (strcmp(s, "--ftr") == 0) {
        print_footer();
      }
      else {
        usage();
      }

      exit(0);
    }
    else if (strcmp(s, "-d") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No depth given\n");
        exit(1);
      }
      else if (depth != 0.0) {
        fprintf(stderr, "Depth given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf%n", &depth, &n) != 1 || o[n] != '\0' || depth <= 0.0) {
        fprintf(stderr, "Wrong depth format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-t") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No times given\n");
        exit(1);
      }
      else if (times != NULL) {
        fprintf(stderr, "Times given twice\n");
        exit(1);
      }

      times = parse_double_list(o, &time_cnt);

      if (times == NULL) {
        fprintf(stderr, "Wrong time format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-g") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No gasses given\n");
        exit(1);
      }
      else if (gasses != NULL) {
        fprintf(stderr, "Gasses given twice\n");
        exit(1);
      }

      gasses = parse_gas_list(o, &gas_cnt);

      if (gasses[0]->is_optional) {
        fprintf(stderr, "The bottom gas should not be optional\n");
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-f") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No gradient factors given\n");
        exit(1);
      }
      else if (gf_min != 0.0) {
        fprintf(stderr, "Gradient factors given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf,%lf%n", &gf_min, &gf_max, &n) != 2 ||
               (o[n] != '\0' && (o[n] != 'r' || o[n + 1] != '\0')) ||
               gf_min <= 0.0 || gf_min > gf_max) {
        fprintf(stderr, "Wrong gradient factor format: %s\n", o);
        exit(1);
      }

      if (o[n] == 'r') {
        smooth_gradient = 0;
      }

      i++;
    }
    else if (strcmp(s, "-s") == 0) {
      int  j;

      if (i + 1 >= argc) {
        fprintf(stderr, "No stops given\n");
        exit(1);
      }
      else if (input_stops != NULL) {
        fprintf(stderr, "Stops given twice\n");
        exit(1);
      }

      input_stops = parse_double_list(o, &input_stop_cnt);

      if (input_stops == NULL) {
        fprintf(stderr, "Wrong stop format: %s\n", o);
        exit(1);
      }

      for (j = 0; j < input_stop_cnt - 2; j++) {
        if (input_stops[j] >= input_stops[j + 1]) {
          fprintf(stderr, "Stops should be increasing in depth: %s\n", o);

          exit(1);
        }
      }

      i++;
    }
    else if (strcmp(s, "-a") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No ambient pressure given\n");
        exit(1);
      }
      else if (p_amb != 0.0) {
        fprintf(stderr, "ambient pressure given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf%n", &p_amb, &n) != 1 || o[n] != '\0' || p_amb <= 0.0) {
        fprintf(stderr, "Wrong ambient pressure format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-r") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No ascent rate given\n");
        exit(1);
      }
      else if (arate != 0.0) {
        fprintf(stderr, "ascent rate given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf%n", &arate, &n) != 1 || o[n] != '\0' || arate <= 0.0) {
        fprintf(stderr, "Wrong ascent rate format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-u") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No gas usage given\n");
        exit(1);
      }
      else if (gas_use_dive != 0.0) {
        fprintf(stderr, "Gas usage given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf,%lf%n", &gas_use_dive, &gas_use_deco, &n) == 2) {
        if (o[n] != '\0' || gas_use_dive <= 0.0 || gas_use_deco <= 0.0) {
          fprintf(stderr, "Wrong gas usage format: %s\n", o);
          exit(1);
        }
      }
      else if (sscanf(o, "%lf%n", &gas_use_dive, &n) != 1 ||
               o[n] != '\0' || gas_use_dive <= 0.0) {
        fprintf(stderr, "Wrong gas usage format: %s\n", o);
        exit(1);
      }

      if (gas_use_deco == 0.0) {
        // only one usage specified
        gas_use_deco  =  gas_use_dive;
      }

      i++;
    }
    else if (strcmp(s, "-p") == 0) {
      if (printer != NULL) {
        fprintf(stderr, "Postscript output specified twice\n");
        exit(1);
      }

      printer = printer_create_ps();
    }
    else if (strcmp(s, "-i") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No initial profile given\n");
        exit(1);
      }
      else if (init_time != 0.0) {
        fprintf(stderr, "Initial profile given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf,%lf%n", &init_depth, &init_time, &n) != 2 || o[n] != '\0' || init_depth < 0.0 || init_time <= 0.0) {
        fprintf(stderr, "Wrong initial profile format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-j") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No finish profile given\n");
        exit(1);
      }
      else if (finish_time != 0.0) {
        fprintf(stderr, "Finish profile given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%lf,%lf%n", &finish_depth, &finish_time, &n) != 2 || o[n] != '\0' || finish_depth < 0.0 || finish_time <= 0.0) {
        fprintf(stderr, "Wrong finish profile format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-k") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No initial gas given\n");
        exit(1);
      }
      else if (init_gas != -1) {
        fprintf(stderr, "Initial gas given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%d%n", &init_gas, &n) != 1 || o[n] != '\0' || init_gas < 0) {
        fprintf(stderr, "Wrong initial gas format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-l") == 0) {
      if (i + 1 >= argc) {
        fprintf(stderr, "No fininsh gas given\n");
        exit(1);
      }
      else if (finish_gas != -1) {
        fprintf(stderr, "Finish gas given twice\n");
        exit(1);
      }
      else if (sscanf(o, "%d%n", &finish_gas, &n) != 1 || o[n] != '\0' || finish_gas < 0) {
        fprintf(stderr, "Wrong finish gas format: %s\n", o);
        exit(1);
      }

      i++;
    }
    else if (strcmp(s, "-z") == 0) {
      if (!total_is_deco) {
        fprintf(stderr, "Full total time set twice\n");
        exit(1);
      }

      total_is_deco = 0;
    }
    else {
      fprintf(stderr, "Unknown option: %s\n", s);
      exit(1);
    }
  }

  if (depth == 0.0) {
    fprintf(stderr, "No depth given\n");
    exit(1);
  }

  if (times == NULL) {
    fprintf(stderr, "No times given\n");
    exit(1);
  }

  if (gasses == NULL) {
    gasses = parse_gas_list("air", &gas_cnt);
  }

  if (gf_min == 0.0) {
    gf_min = GF_MIN;
    gf_max = GF_MAX;
  }

  if (p_amb == 0.0) {
    p_amb = P_STD;
  }

  if (arate == 0.0) {
    arate = 10.0;
  }

  if (printer == NULL) {
    printer = printer_create_human();
  }

  if (input_stop_cnt >= 2 && input_stops[input_stop_cnt - 2] > input_stops[input_stop_cnt - 1]) {
    stop_step  =  input_stops[input_stop_cnt - 1];
    input_stop_cnt--;
  }

  schedule  =  schedule_create(input_stop_cnt, gas_cnt, p_amb, gf_min, gf_max, smooth_gradient, arate);

  schedule->gas_use_dive  =  gas_use_dive;
  schedule->gas_use_deco  =  gas_use_deco;

  for (i = 0; i < input_stop_cnt; i++) {
    schedule->fixed_stops[i] = input_stops[i];
  }

  free(input_stops);

  schedule->stop_step  =  stop_step;

  for (i = 0; i < gas_cnt; i++) {
    schedule->gasses[i]  =  gasses[i];
  }

  if (init_gas == -1) {
    init_gas = 0;
  }

  if (init_gas >= gas_cnt) {
    fprintf(stderr, "Initial gas index too high: %d\n", init_gas);
    exit(1);
  }

  if (finish_gas >= gas_cnt) {
    fprintf(stderr, "Finish gas index too high: %d\n", finish_gas);
    exit(1);
  }

  schedule->init_time     =  init_time;
  schedule->init_depth    =  init_depth;
  schedule->init_gas      =  init_gas;

  schedule->finish_time   =  finish_time;
  schedule->finish_depth  =  finish_depth;
  schedule->finish_gas    =  finish_gas;

  free(gasses);

  context                 =  (context_t*) malloc(sizeof(context_t));
  context->schedule       =  schedule;
  context->depth          =  depth;
  context->time_cnt       =  time_cnt;
  context->times          =  times;
  context->printer        =  printer;
  context->total_is_deco  =  total_is_deco;

  return  context;
}


void  context_free( context_t*  context )
{
  printer_free(context->printer);
  free(context->times);
  schedule_free(context->schedule);
  free(context);
}


table_t*  table_create( par_t*       par,
                        schedule_t*  schedule,
                        double       depth,
                        double*      times,
                        int          time_cnt )
{
  state_t*  state       =  state_create(schedule, par);
  state_t*  tmp_state   =  state_create(schedule, par);
  table_t*  table       =  (table_t*) malloc(sizeof(table_t));
  gas_t*    air         =  gas_parse("air", 0, NULL);
  double    deco_depth  =  depth;
  double    p;
  int       i;

  if (schedule->finish_time != 0.0 && schedule->finish_depth > depth) {
    deco_depth  =  schedule->finish_depth;
  }

  table->time_cnt   =  time_cnt;
  table->stop_list  =  stop_list_create(schedule, deco_depth);
  table->decos      =  (deco_t**) malloc(time_cnt * sizeof(deco_t*));

  for (i = 0; i < time_cnt; i++) {
    int  j;

    saturate(schedule, par, state, schedule->p_amb, air);

    if (schedule->init_time != 0.0) {
      p = depth_to_p(schedule, schedule->init_depth);

      dive(schedule, par, state, p, p, schedule->init_time, schedule->init_gas);
    }

    p = depth_to_p(schedule, depth);

    dive(schedule, par, state, p, p, times[i], 0);

    table->decos[i] = get_deco(par, state, schedule, deco_depth, table->stop_list, tmp_state);
    table->decos[i]->total_time = state->time;

    for (j = 0; j < schedule->gas_count; j++) {
      table->decos[i]->gas_use[j] = state->gas_use[j];
    }
  }

  gas_free(air);
  state_free(state);
  state_free(tmp_state);

  return  table;
}


void  table_free( table_t*  table )
{
  int  i;

  for (i = 0; i < table->time_cnt; i++) {
    deco_free(table->decos[i]);
  }

  stop_list_free(table->stop_list);
  free(table->decos);
  free(table);
}


int  main( int     argc,
           char**  argv )
{
  context_t*    context    =  parse_options(argc, argv);
  schedule_t*   schedule   =  context->schedule;
  double        depth      =  context->depth;
  int           time_cnt   =  context->time_cnt;
  par_t*        par        =  par_create();
  table_t*      table0;

  context->printer->start_table(time_cnt + 2);

  table0 = table_create(par, schedule, depth, context->times, context->time_cnt);
  table_print(context->printer, schedule, table0->stop_list, NULL, depth, context->times, table0->decos, time_cnt, context->total_is_deco);

  int i;

  for (i = schedule->gas_count - 1; i >= 1; i--) {
    schedule_t*  new_schedule  =  schedule_derive(schedule, i);
    table_t*     table;

    if (new_schedule == NULL) {
      continue;
    }

    context->printer->drule();

    table = table_create(par, new_schedule, depth, context->times, context->time_cnt);
    table_print(context->printer, new_schedule, table->stop_list, table0->stop_list, depth, context->times, table->decos, time_cnt, context->total_is_deco);
    table_free(table);

    schedule_free(new_schedule);
  }

  context->printer->end_table();

  table_free(table0);
  par_free(par);
  context_free(context);

  return  0;
}
