#include "clucu_error.h"
#include <string.h> //strrchr()函数所需头文件
//linux :
#define filename(x) strrchr(x,'/')?strrchr(x,'/')+1:x


// Convenience function to handle warnings

void clucu_raise_warning1(int err,const char *file,const char *function, int line,const char* msg, ...) 
{
    char buffer[256];

    va_list va;
    va_start(va, msg);//取第一个 可变参数 的指针给va, msg是函数声明中的最后一个 固定参数。
    vsnprintf(buffer, 250, msg, va);
    //将格式化数据从 可变参数列表va 写入 缓冲区buffer
    //vsnprintf会自动在写入字符的后面加上停止符\0。如上，缓存区buffer最大256个字符，在写入250个字符后，自动添加了\0。
    va_end(va);
   fprintf(stderr, "WARNING %d: %s/%s()%d: %s\n", err,filename(file),function,line, buffer);
   abort();//abort函数用于不正常地终止一个正在执行的程序.它可能不会清理包含未输出数据的输出缓冲区，不会关闭打开的流，也不会删除临时文件
   //exit(1);//exit()函数用于正常终止程序。exit(0) 表示程序正常退出,exit⑴/exit(-1)表示程序异常退出。
}
void clucu_raise_warning(int err,const char* msg, ...) {
    char buffer[256];

    va_list va;
    va_start(va, msg);//取第一个 可变参数 的指针给va, msg是函数声明中的最后一个 固定参数。
    vsnprintf(buffer, 250, msg, va);
    //将格式化数据从 可变参数列表va 写入 缓冲区buffer
    //vsnprintf会自动在写入字符的后面加上停止符\0。如上，缓存区buffer最大256个字符，在写入250个字符后，自动添加了\0。
    va_end(va);

    // For now just print warning to stderr if debug is enabled.
    // TODO: Implement some kind of error stack that can be passed on to, e.g.,
    // the python binding.
    /*
    if (_clucu_debug_mode_policy == CLUCU_DEBUG_MODE_ON) {
        
    }
    */
   fprintf(stderr, "WARNING %d: %s\n", err, buffer);
   //abort();//abort函数用于不正常地终止一个正在执行的程序.它可能不会清理包含未输出数据的输出缓冲区，不会关闭打开的流，也不会删除临时文件
   //exit(1);//exit()函数用于正常终止程序。exit(0) 表示程序正常退出,exit⑴/exit(-1)表示程序异常退出。
}
// Convenience function to handle warnings
void clucu_raise_gsl_warning(int gslstatus, const char* msg, ...) 
{
    char buffer[256];

    va_list va;
    va_start(va, msg);
    vsnprintf(buffer, 250, msg, va);
    va_end(va);

    clucu_raise_warning(gslstatus, "%s: GSL ERROR: %s\n", buffer, gsl_strerror(gslstatus));
    //exit(1);//exit()函数用于正常终止程序。
}
void clucu_raise_gsl_warning1(int gslstatus,const char *file,const char *function, int line,const char* msg, ...)
{
    char buffer[256];

    va_list va;
    va_start(va, msg);
    vsnprintf(buffer, 250, msg, va);
    va_end(va);
    fprintf(stderr, "WARNING %d: %s/%s()%d: %s GSL ERROR: %s:\n", gslstatus,filename(file),function,line, buffer,gsl_strerror(gslstatus));

    //clucu_raise_warning1(gslstatus,file,function,line, "%s: GSL ERROR: %s\n", buffer, gsl_strerror(gslstatus));
    //exit(1);//exit()函数用于正常终止程序。
}