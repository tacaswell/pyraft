#ifndef RAFT_WRAPPER_XFCTFROMSPEC_H
#define RAFT_WRAPPER_XFCTFROMSPEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../raft/raft_image.h"

void getFluorFromSpec(raft_image sino,
		      raft_image table, //spec table from XFCT 
	              int s,      //slice number 
		      int x,      //number of points: x 
		      int y,      //number of points: y
		      int a,      //number of angles: a
		      int t,      //total number of columns @ file
		      int col);   //column to extract fluorescence data

void getTransFromSpec(raft_image sino,
		      raft_image table, //spec table from XFCT 
	              int s,      //slice number 
		      int x,      //number of points: x 
		      int y,      //number of points: y
		      int a,      //number of angles: a
		      int t,      //total number of columns @ file
		      int c,	  //column to extract transmission data
		      int n);	  //column to extract normalization factor

#ifdef __cplusplus
} //extern "C" {
#endif

#endif // #ifndef RAFT_WRAPPER_XFCTFROMSPEC_H

