#ifndef _OPED_H_
#define _OPED_H_

#include "interp.h"
#include <math.h>
#include "raft_param.h"
#include "raft_oped.h"
#include "raft_scan.h"


double 
eval_opedsum(raft_scan_t *data,
	     raft_oped_t *workspace,
	     int j,
	     int k);



#endif
