# LocalitySensitiveHashing
Execution Instructions:

1) Move the file to the target location
2) Run the command “mpicc parallelKmeans.c -o parallelKmeans -lm” to compile
3) Execute with the command “mpirun -np noofprocesses ./parallelKmeans arg1 arg2 arg3 ”
4) Input number of processes as “noofprocesses”
3) Input number of Dimensions as “arg1”
4) Input number of Points as “arg2”
5) Input number of Clusters as “arg3”
6) Output is the number of visits to other points and nearest neighbor distance, for 10 search points.
