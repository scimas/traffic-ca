all : traffic_serial.c traffic_parallel.c
	gcc traffic_serial.c -o e_traffic_serial -Wall -O3 -lm
	gcc traffic_parallel.c -o e_traffic_parallel -Wall -O3 -lm -fopenmp

serial : traffic_serial.c
	gcc traffic_serial.c -o e_traffic_serial -Wall -O3 -lm

parallel : traffic_parallel.c
	gcc traffic_parallel.c -o e_traffic_parallel -Wall -O3 -lm -fopenmp

clean :
	rm e_traffic_parallel e_traffic_serial
