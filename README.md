
## setup 
- install the python env provided in the yml file
- the environment variable `$PATH_SRC_PYTHON_LOCAL` need to be defined and point where you have the dir `Process_MNH`
```
mamba env create -f environment.yml
```

## run the script
```
source ~/miniforge3/bin/activate process_mnh
python do_processMNH.py -i PATH_WITH_NC_FILES -t 0 -f nc
```
for more details on the input
```
python do_processMNH.py  --help
```

if you want to add lambda2 criteria in the output run :
```
python do_processMNH.py -i /data/paugam/MNH/IdealCase3/004_ffb-mesonh/ -var lambda2,lambda2
```


## old stuff
To compile writeplt on kraken
```
module pugre
module load compiler/gcc/5.4
module load lib/netcdf-c/4.6.1_gcc
make python
```
last command failed
run it manually with the `-shared` option
```
/home/logiciels/local/gcc/5.4.0/bin/gfortran -shared -Wall -g /home/globc/paugam/Src/Process_MNH/src/writePlt/lib/libtecio.a -lstdc++ -L/softs/local/netcdf-c/4.6.1/lib -lnetcdf ./src.linux-x86_64-2.7/_writePltmodule.o ./src.linux-x86_64-2.7/src.linux-x86_64-2.7/fortranobject.o ./f90wrap_cons.o ./f90wrap_save_tecplot.o ./f90wrap_toplevel.o lambda2.o save_tecplot.o tred1.o r1mach.o pythag.o cons.o tqlrat.o /home/globc/paugam/Src/Process_MNH/src/writePlt/lib/libtecio.a -L/softs/local/netcdf-c/4.6.1/lib -L/home/logiciels/local/gcc/5.4.0/bin/../lib/gcc/x86_64-unknown-linux-gnu/5.4.0/../../../../lib64 -L/home/logiciels/local/gcc/5.4.0/bin/../lib/gcc/x86_64-unknown-linux-gnu/5.4.0/../../../../lib64 -L/home/globc/paugam/anaconda2/envs/mypy_cpu/lib -lstdc++ -lnetcdf -lpython2.7 -lgfortran -o ./_writePlt.so
```
