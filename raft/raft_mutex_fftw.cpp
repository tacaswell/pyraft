#include "raft_mutex_fftw.h"
#include <mutex>

std::mutex raft_fftw_mutex;

void raft_mutex_fftw_lock( void )
{
   raft_fftw_mutex.lock();
}

void raft_mutex_fftw_unlock( void )
{
   raft_fftw_mutex.unlock();
}
