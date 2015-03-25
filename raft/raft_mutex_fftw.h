#ifndef RAFT_MUTEX_FFTW_H
#define RAFT_MUTEX_FFTW_H

#ifdef __cplusplus
extern "C" {
#endif

// These are function used to protect plan creation
// routines from fftw to misbehave in parallel usage.

void raft_mutex_fftw_lock( void );
void raft_mutex_fftw_unlock( void );

#ifdef __cplusplus
} // extern "C"
#endif

#endif // #ifndef RAFT_FFTW_MUTEX_H
