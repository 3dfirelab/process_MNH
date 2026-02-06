#/bin/bash
make clean
f2py -m interp -h interp.pyf src/*.f90 --overwrite-signature
f2py -c --fcompiler=gnu95 interp.pyf src/*.f90  --f90flags='-fdefault-real-8 -fdefault-double-8 -fcray-pointer -O0 -g -fbounds-check -Wall -fbacktrace'
