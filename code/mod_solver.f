#include "solver_comm.h"

	module mod_solver
	implicit none
	save        

	!======================================================================

	real(8), parameter :: tol = 1.0E-5

	integer :: nstep_max     = 300000
	integer :: step4out_u    = 20000
	integer :: step4out_iter = 10000

	integer :: iter 
	integer :: len_recl 

	real(8) :: omega
	real    :: time_solver 

!======================================================================
		
      end module mod_solver