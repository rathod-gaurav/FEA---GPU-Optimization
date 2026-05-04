## Running the Solver
Mesh resolution and thread count are set via compile-time constants (or environment variables via run.sh). The benchmark mesh sizes (20×6×6, 40×12×12, 80×24×24, 160×48×48) can be selected by editing the mesh parameters in main.cpp or passing them as arguments (see run.sh). Output VTK files are written to the output/ directory and can be visualized in ParaView.
The following steps build and run the optimized SVK/NH/MR solvers:
cd SVK\ Solver\ -\ Parallel/
sbatch run.sh
- - - - - - - - - 
cd NH\ Solver\ -\ Parallel/
sbatch run.sh
- - - - - - - - - 
cd MR\ Solver\ -\ Parallel/
sbatch run.sh

To reproduce the scaling results reported in Sections 4.2 and 4.3:
Change directory to one of the above solvers’ folders.
sbatch experiment20.sh
sbatch experiment40.sh
sbatch experiment80.sh
sbatch experiment160.sh
The output files will be stored inside experiments/experiment{20/40/80/160} folders inside respective solvers’ folders. The last line of each output file will show the total time taken by the program to execute.
