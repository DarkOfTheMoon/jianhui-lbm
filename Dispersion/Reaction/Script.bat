@echo off

set CodeFolder=E:\Reaction\Codes
set InputFolder=E:\Reaction\Input
set OutputFolder=E:\Reaction\Output
set WorkingFolder=E:\Reaction\Working
set GeometryFolder=E:\Reaction\Geometry

set LBExe="D:\Visual Studio 2010 Projects\Lattice Boltzmann\x64\Release\Lattice Boltzmann.exe"
set ReactionExe="D:\Visual Studio 2010 Projects\Reaction\x64\Release\Reaction.exe"

set LBInputFile="E:\Reaction\Input_LB_RW.txt"

set LBThreads=8


del %WorkingFolder%\* /Q


copy %LBExe% %CodeFolder%\LB.exe /Y
copy %ReactionExe% %CodeFolder%\Reaction.exe /Y

copy %GeometryFolder%\geo.bin %InputFolder%\geo.bin /Y
copy %GeometryFolder%\vel.bin %InputFolder%\vel.bin /Y

copy %GeometryFolder%\geo.bin %WorkingFolder%\geo.bin /Y
copy %GeometryFolder%\vel.bin %WorkingFolder%\vel.bin /Y


FOR /L %%i IN (0 1 200) DO (


%CodeFolder%\Reaction.exe


copy %WorkingFolder%\geo.bin %CodeFolder%\geo.bin /Y
copy %WorkingFolder%\vel.bin %CodeFolder%\vel.bin /Y
copy %WorkingFolder%\GeometryUpdate.txt GeometryUpdate.txt /Y


mpiexec -np %LBThreads% %CodeFolder%\LB.exe %LBInputFile% 1


copy %CodeFolder%\geo.bin %WorkingFolder%\geo.bin /Y
copy %CodeFolder%\vel.bin %WorkingFolder%\vel.bin /Y

)

:END

PAUSE