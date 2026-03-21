# g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -Wall -O3 -std=c++17 -lcurl -o main
# # /mnt/c/Users/ratho/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp 
# ./main


# gprof
g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -pg -O2 -lcurl -o main_gprof
./main_gprof
gprof main_gprof gmon.out > analysis_gprof.txt

# perf
g++ -I /home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0 main.cpp -g -O2 -lcurl -o main_perf
perf record -g ./main_perf
perf report #visual report in the terminal 
perf stat ./main_perf #a top-down report of the program's performance, including metrics such as CPU cycles, instructions, cache misses, etc.
perf report --stdio > perf_report.txt #save the report to a text file