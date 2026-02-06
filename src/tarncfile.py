import subprocess
import glob
import os
import sys
import shutil

currentdir = os.getcwd()
tmp = currentdir.split('/')
tarfilename = '{:s}_{:s}.tar.gz'.format(tmp[-2],tmp[-1])

ncfiles = sorted(glob.glob('Postproc/*'))

listargs = ['tar', '--overwrite', '-cvzf', tarfilename,]
for ncfile_ in ncfiles:
    listargs.append(ncfile_)
subprocess.call(listargs)

#shutil.move(tarfilename, './Postproc/outputFile/Netcdf/')
