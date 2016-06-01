#ifndef _RAFT_ERRNO_H_
#define _RAFT_ERRNO_H_

/*######################################################
  Title: Error number
  ####################################################*/

/*+====================================================+
  
  CONSTANT: enum

  RAFT error numbers.

  RAFT_SUCCESS - successfull operation
  RAFT_ENOMEM - not enough memory
  RAFT_EDOM - domain error
  RAFT_CFG - description file error
  
  +====================================================+
*/

enum{
  RAFT_SUCCESS = 1000,
  RAFT_ENOMEM = 1001,
  RAFT_EDOM = 1002,
  RAFT_CFG = 1003
};

#endif
