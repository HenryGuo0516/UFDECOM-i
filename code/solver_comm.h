!
! preprocess
!
! compile with option -fpp -cpp
!
#define int integer(4)
#define long integer(8)
#define float real(4)
#define cdouble real(4)
#define char character
#define bool logical

! ======================================
! parameters for testing matrix
! #define _N8  ! matrix of 8-by-8
#ifdef _N8

#define N 8
#define nnz 22 ! (3*N-2) 
#define _C_DATA  ! read data from C/C++
!#define itermax N

#else

#define N 43254
#define nnz 130816
!#define itermax 1000

#endif
!====================================== 

