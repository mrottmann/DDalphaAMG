/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#ifndef THREADING_H
#define THREADING_H


#ifdef FLAT_OMP

// flat omp: do not distinguish between threads on cores and hyperthreads
#define CORE_BARRIER(threading) \
    do { \
    threading->thread_barrier(threading->thread_barrier_data, threading->core/60); \
    if(threading->core/60 == 0) \
        threading->barrier(threading->core); \
    threading->thread_barrier(threading->thread_barrier_data, threading->core/60); \
    } while(0)
#define HYPERTHREAD_BARRIER(threading) \
    do { \
    threading->thread_barrier(threading->thread_barrier_data, threading->core/60); \
    } while(0)

#else

// nested omp: split into cores, each core splits into hyperthreads (like DD preconditioner)
#define CORE_BARRIER(threading) \
    do { \
        threading->barrier(threading->core); \
    } while(0)
#define HYPERTHREAD_BARRIER(threading) \
    do { \
    threading->thread_barrier(threading->thread_barrier_data, threading->thread); \
    } while(0)

#endif


#define START_UNTHREADED_FUNCTION(threading) \
    if(threading->thread != 0) \
        return; \
    CORE_BARRIER(threading); \
    if(threading->core != 0) \
    { \
        CORE_BARRIER(threading); \
        return; \
    }
#define END_UNTHREADED_FUNCTION(threading) \
    CORE_BARRIER(threading);

// only one thread (master) will execute the code section between
// START_LOCKED_MASTER and END_LOCKED_MASTER, and it is protected by barriers
// among cores to prevent data races
#define START_LOCKED_MASTER(threading) \
    if(threading->thread == 0) \
        CORE_BARRIER(threading); \
    if(threading->core + threading->thread == 0) {
#define END_LOCKED_MASTER(threading) \
    } \
    if(threading->thread == 0) \
        CORE_BARRIER(threading);

#define START_MASTER(threading) \
    if(threading->core + threading->thread == 0) {
#define END_MASTER(threading) \
    }

#define SYNC_MASTER_TO_ALL(threading) \
    if(threading->thread == 0) \
        CORE_BARRIER(threading); \
    HYPERTHREAD_BARRIER(threading);
    
#define SYNC_CORES(threading) \
    if(threading->thread == 0) \
        CORE_BARRIER(threading);
#define SYNC_HYPERTHREADS(threading) \
    HYPERTHREAD_BARRIER(threading);

#define START_NO_HYPERTHREADS(threading) \
    if(threading->thread == 0) {
#define END_NO_HYPERTHREADS(threading) \
    }


#ifdef OPENMP
#include <omp.h>
#else
static inline int omp_get_thread_num( void ) {
  return 0;
}
static inline int omp_get_num_threads( void ) {
  return 1;
}
#endif

struct level_struct;

struct common_thread_data
{
    // barrier among cores
    void (*barrier)(int);
    // barrier among hyperthreads on a core
    void (*thread_barrier)(void *, int);
    // *common* workspace for *all* threads
    // sometimes threads need to exchange data, they can use this
    char *workspace;
};

void init_common_thread_data(struct common_thread_data *common);

// holds information relevant for specific core/thread
typedef struct Thread
{
    int core;
    int n_core;
    // for SMT/hyperthreading: threads per core (1-4 on KNC)
    int thread;
    int n_thread;

    // level_struct.num_inner_lattice_sites is split among cores
    // these variables define start and end site for this specific *core* (not thread)
    // but num_inner_lattice_sites depends on the level
    // use level_struct.depth as index
    int start_site[4];
    int end_site[4];
    int n_site[4];
    // index = site*num_lattice_site_var = inner_vector_size
    int start_index[4];
    int end_index[4];
    int n_index[4];

    // barrier among cores
    void (*barrier)(int);
    // barrier among hyperthreads on a core
    void (*thread_barrier)(void *, int);
    void *thread_barrier_data;

    // *common* workspace for *all* threads
    // sometimes threads need to exchange data, they can use this
    char *workspace;
} Thread;

void setup_threading(struct Thread *threading, struct common_thread_data *common, struct level_struct *l);
/* external means the caller gives us all info about threads, and is responsible to later set a proper barrier */
void setup_threading_external(struct Thread *threading, struct common_thread_data *common, struct level_struct *l,
        int n_core, int n_thread, int core, int thread);
void update_threading(struct Thread *threading, struct level_struct *l);
void setup_no_threading(struct Thread *no_threading, struct level_struct *l);

// computes start and end indices for a core inside an array
// puts zero for other hyperthreads
void compute_core_start_end(int start, int end, int *core_start, int *core_end,
        struct level_struct *l, struct Thread *threading);
void compute_core_start_end_custom(int start, int end, int *core_start, int *core_end,
        struct level_struct *l, struct Thread *threading, int granularity);

void finalize_common_thread_data( struct common_thread_data *common );

void finalize_no_threading( struct Thread *no_threading );

extern struct Thread *no_threading;

extern int threaded;

#endif // THREADING_H
