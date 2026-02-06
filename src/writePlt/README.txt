on kraken
when doign make the shared option is missing 
you have to run it manually 

make -n python to show command that you need to run

/usr/bin/gfortran -shared -Wall -g -Wl,-rpath,/lib -L/lib ./src.linux-x86_64-2.7/_writePltmodule.o ./src.linux-x86_64-2.7/src.linux-x86_64-2.7/fortranobject.o ./f90wrap_cons.o ./f90wrap_save_tecplot.o ./f90wrap_toplevel.o lambda2.o save_tecplot.o tred1.o r1mach.o pythag.o cons.o tqlrat.o /home/globc/paugam/Src/Process_MNH/src/writePlt/lib/libtecio.a -L/softs/local/netcdf-c/4.6.1/lib -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/home/globc/paugam/anaconda2/envs/mypy_cpu/lib -lstdc++ -lnetcdf -lpython2.7 -lgfortran -o ./_writePlt.so
