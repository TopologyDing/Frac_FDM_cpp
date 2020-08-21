#pragma once
#include <iostream>

#ifndef NDEBUG
#   define MyAssert(Expr, Msg)  _M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define MyAssert(Expr, Msg) ;
#endif

void _M_Assert(const char* expr_str, bool expr, const char* file, int line, const char* msg)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t" << msg << "\n"
            << "Expected:\t" << expr_str << "\n"
            << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}