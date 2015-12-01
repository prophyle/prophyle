# distutils: language = c++

#from libcpp.unordered_map cimport unordered_map
#from libcpp.utility cimport pair

cdef extern from "cfib.h":
    int _fib "fib"(char *fasta)

def fib(n):
    ''' Returns the nth Fibonacci number.'''
    return _fib("a.fa")
