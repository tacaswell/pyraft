#ifndef _SCAN_H_
#define _SCAN_H_

#include "raft_param.h"
#include <string.h>
#include <confuse.h>


int
parse_geometry(cfg_t *section,
	       int *nrays,
	       int *nviews,
	       int *size,
	       double *min,
	       double *max,
	       int *scantype);


int 
parse_energy(cfg_t *section,
	     raft_scan_t *data);
  


int 
parse_statistic(cfg_t *section,
		raft_scan_t *data);



int
parse_mesh(cfg_t *section,
	   raft_scan_t *data);


#endif
