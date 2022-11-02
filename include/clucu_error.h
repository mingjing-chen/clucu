#pragma once
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
//这里，如果可变参数被忽略或为空，'##'操作将使预处理器(preprocessor)去除掉它前面的那个逗号。
#define CLUCU_RAISE_WARNING(err,...) \
       clucu_raise_warning1 (err, __FILE__,__FUNCTION__, __LINE__,##__VA_ARGS__) ; 

#define CLUCU_RAISE_GSL_WARNING(gslstatus,...) \
       clucu_raise_gsl_warning1 (gslstatus, __FILE__,__FUNCTION__, __LINE__,##__VA_ARGS__) ; 

#define CLUCU_CHECK_ERROR(x) if(*status){CLUCU_RAISE_WARNING(*status,NULL);return x;}

//malloc, gsl的alloc
#define CLUCU_ERROR_MEMORY 1025
#define CLUCU_ERROR_LINSPACE 1026
#define CLUCU_ERROR_LOGSPACE 1052
#define CLUCU_ERROR_LOG10SPACE 24231244
#define CLUCU_ERROR_LINLOGSPACE 1053
#define CLUCU_ERROR_FILE 23453453
#define CLUCU_ERROR_COMPUTE 35235232
#define CLUCU_ERROR_INCONSISTENT 1027
#define CLUCU_ERROR_SPLINE 1028
#define CLUCU_ERROR_SPLINE_EV 1029
#define CLUCU_ERROR_ROOT 1031
#define CLUCU_ERROR_INTEG 1030
#define CLUCU_ERROR_COMPUTECHI 1033
#define CLUCU_ERROR_PARAMETERS 1036
#define CLUCU_ERROR_NU_INT 1037
#define CLUCU_ERROR_NOT_IMPLEMENTED 1040
#define CLUCU_ERROR_MNU_UNPHYSICAL 1041


#define CLUCU_ERROR_GROWTH_INIT 1058
#define CLUCU_ERROR_DISTANCES_INIT 1059
#define CLUCU_ERROR_SIGMA_INIT 1062
#define CLUCU_ERROR_OVERWRITE 1064

void clucu_raise_warning1(int err,const char *file,const char *function, int line,const char* msg, ...) ;

/** Raise a warning
 * Given a status, give a warning message.
 * @return void
 */
void clucu_raise_warning(int err, const char* msg, ...);

/** Raise a warning based on a GSL error message
 * Given a GSL status, give a warning message.
 * @return void
 */
void clucu_raise_gsl_warning(int gslstatus, const char* msg, ...);


void clucu_raise_gsl_warning1(int gslstatus,const char *file,const char *function, int line,const char* msg, ...);
