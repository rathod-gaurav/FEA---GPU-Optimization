# g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -Wall -O3 -std=c++17 -lcurl -o main
# # /mnt/c/Users/ratho/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp 
# ./main


# gprof
g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -pg -O2 -lcurl -o main_gprof
./main_gprof
gprof main_gprof gmon.out > analysis_gprof.txt

# perf
g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -g -O2 -lcurl -o main_perf
perf record -g --all-user taskset -c 3 ./main_perf #record the performance data of the program, including call graphs, while running it on CPU core 3. The --all-user flag ensures that only user-space code is profiled, excluding kernel code.
perf report #visual report in the terminal 
perf stat ./main_perf #It will rerun the code. a top-down report of the program's performance, including metrics such as CPU cycles, instructions, cache misses, etc.
perf report --stdio > perf_report_isolated_core.txt #save the report to a text file

#intel vtune
# source /opt/intel/oneapi/setvars.sh
# Always compile with -g (debug symbols) + optimizations kept on
g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -O2 -g -lcurl -march=native -o main_vtune
# Run the VTune profiler with the appropriate analysis type (e.g., hotspots, memory access, etc.)
vtune -collect performance-snapshot hotspots memory-access hpc-performance \
      -knob collect-memory-bandwidth=true \
      -r ./roofline_result \
      -- ./main_vtune

vtune -report summary -r ./roofline_result --format text > vtune_report.txt #Generate a summary report in text format and save it to vtune_report.txt
vtune-gui ./roofline_result #Open the VTune GUI to visualize the results in roofline_result directory