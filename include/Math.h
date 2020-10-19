/*********************************************************************************
 * 
 *                  Math.h (Morfo70)
 * 
 *  auxiliar functions to solve tridiagonal systems A*y=f (periodic o no)
 *  Implementation in:
 *      - Math.cpp 
 * 
 * *******************************************************************************/
#ifndef MATH_H_
#define MATH_H_
#include <iostream>

namespace Math
{   
    /**   Solve the  n x n  tridiagonal system for y:
    *   [ a(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
    *   [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
    *   [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
    *   [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
    *   [                    ...    ...    ...        ] [        ]   [        ]
    *   [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
    *   [                                 b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
    *
    *   f must be a vector (row or column) of length n
    *   a, b, c must be vectors of length n (note that b(1) and c(n) are not used)
    *   This funtion solves M linear systems simulstaneally, length(Â)=length(B)....=n*m */
    void TridiagonalSolver(double  * Y,double  * V, double* A,  double* B,  double* C, double * F, int n,int m);
    
    
    /**   Solve a periodic tridiagonal system using El-Mikkawy 2005 algorithm
        *   Solve the  n x n  tridiagonal system for y:
        *   
        *   [ a(1)  c(1)                            b(1)  ] [  y(1)  ]   [  f(1)  ]
        *   [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
        *   [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
        *   [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
        *   [                    ...    ...    ...        ] [        ]   [        ]
        *   [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
        *   [c(n)                             b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]

        *   f must be a vector (row or column) of length n
        *   a, b, c must be vectors of length n */
    void TridiagonalSolverPeriodic(double* a,  double* b,  double* c, double * f,  int n,double * y);
    
    /**   Solve a periodic tridiagonal system using El-Mikkawy 2005 algorithm
    *   Solve the  n x n  tridiagonal system for y:
    * 
    *   [ a(1)  c(1)                            b(1)  ] [  y(1)  ]   [  f(1)  ]
    *   [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
    *   [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
    *   [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
    *   [                    ...    ...    ...        ] [        ]   [        ]
    *   [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
    *   [c(n)                             b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
    * 
    *   f must be a vector (row or column) of length n
    *   a, b, c must be vectors of length n 
    *   This funtion solves M linear systems simulstaneally, length(Â)=length(B)....=n*m */
    void TridiagonalSolverPeriodic(double * Y, double* A,  double* B,  double* C, double * F, int n, int m);  
    
 }
#endif