# Performance Test #

In this folder is contained a script to best run a performance test of the `solver` class. <br>

It compiles and links the `performance_test` executable in the `../src/` folder, then runs it with
1, 2 and 4 processes. `performance_test` solves the default problem shown in the other `README.md`
file on meshes using 16, 32, 64, 128 and 256 points per side. <br>

For each mesh refinement level, the solution is calculated `REPS` times to get the average amount 
of time spent calculating it, where `REPS` is a pre-processor variable set to 1 by default. <br>
The results and statistics of each test are stored in the data folder. <br>

`plot_results.py` is then called by the script to create plots of the time taken and $L^2$ error
as the mesh refinement increases. <br>

The script is run from the command line as `./performance_test.sh REPS`. The value of `REPS` in the last
call of the script is stored inside `cache.dat` so that the `performance_test` executable only gets
rebuilt if the new call of the script has a different value of `REPS`.<br>


## Warnings ##

Using the `Makefile` in the `../src/` folder with a manually defined `REPS` (such as
`make performance_test -DREPS=4`) will **not** update `cache.dat`. Using a simple `make performance_test` 
or just `make` will, on the other hand, read the value of `REPS` from `cache.dat`, so that no problems arise. <br>

Another possible problem would arise from updating one of the dependencies of the `performance_test` executable 
**and** using the same value of `REPS` as the one in `cache.dat` without updating or removing the executable first, as this
will **not** rebuild it.