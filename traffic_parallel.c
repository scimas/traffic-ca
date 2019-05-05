#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

typedef struct cars {
    int x;                  // x position on the road
    int y;                  // Lane number on the road
    int v;                  // Current velocity
    int v_d;                // Desired velocity
    int lane_change_now;    // Did lane change this iteration
    int lane_change_prev;   // Did lane change last iteration
} car;

typedef struct {
    int LANES;          // Number of lanes (currently only 2 possible)
    int L;              // Length of the road
    int N;              // Number of cars
    int MAX_ITER;       // Iterations
    int v_max;          // Maximum allowed velocity
    double p_change;    // Lane change probability
    double p_break;     // Random braking probability
} parameters;

void initialise(car *cars, int *grid, int *new_grid, parameters sim, int seed);
void update(car *cars, int *grid, int *new_grid, parameters sim, int seed, int iter);
car * preceedingCar(car *cars, int car_index, int *grid, parameters sim);
int aheadThisLane(car *cars, int car_index, parameters sim);
int aheadOtherLane(car *cars, int car_index, parameters sim);
int behindOtherLane(car *cars, int car_index, parameters sim);
int gapAhead(car *cars, int car_index, int *grid, parameters sim);
int gapAheadOther(car *cars, int car_index, int *grid, parameters sim);
int gapBehindOther(car *cars, int car_index, int *grid, parameters sim);
void printGrid(car *cars, int *grid, int length, parameters sim);
int mod(int a, int n);

int main(int argc, char **argv) {
    if (argc == 1) {
        printf("Usage: ./traffic_serial seed\n");
        printf("seed: integer used to seed the PRNG\n");
        exit(99);
    }

    parameters sim;
    sim.LANES = 2;
    sim.L = 133333;
    sim.N = (int) sim.LANES * sim.L * 0.08;
    sim.MAX_ITER = 5000;
    sim.v_max = 5;
    sim.p_change = 1.0;
    sim.p_break = 0.1;

    car *cars = malloc(sizeof *cars * sim.N);
    int *grid = malloc(sizeof *grid * sim.LANES * sim.L);
    int *new_grid = malloc(sizeof *new_grid * sim.LANES * sim.L);

    if (cars == NULL || grid == NULL || new_grid == NULL) {
        printf("Could not allocate memory!!!\n");
        exit(1);
    }

    double density = (double) sim.N / (double) (sim.LANES * sim.L);
    double flow = 0;
    double lane_changes = 0;
    double ping_pong_changes = 0;

    int seed = atoi(argv[1]);

    initialise(cars, grid, new_grid, sim, seed);
    printf("Initialisation complete.\n\n\n");

    // struct timespec sleeper = {0, 400000000};
    // struct timespec remain;

    for (int i = 0; i < 1000; i++) {
        update(cars, grid, new_grid, sim, seed, i);
        // printGrid(cars, grid, 100, sim);
        // nanosleep(&sleeper, &remain);
    }

    for (int i = 0; i < sim.MAX_ITER; i++) {
        update(cars, grid, new_grid, sim, seed, i);
        #pragma omp parallel for reduction(+:lane_changes, ping_pong_changes)
        for (int j = 0; j < sim.N; j++) {
            lane_changes += cars[j].lane_change_now;
            if (cars[j].lane_change_now * cars[j].lane_change_prev == 1) {
                ping_pong_changes += 1;
            }
        }
        if (i % 5 == 0) {
            #pragma omp parallel for reduction(+:flow)
            for (int j = 0; j < sim.N; j++) {
                flow += cars[j].v;
            }
        }
    }
    flow /= (double) (sim.LANES * sim.L * sim.MAX_ITER);
    lane_changes /= (double) (sim.LANES * sim.L * sim.MAX_ITER);
    ping_pong_changes /= (double) (sim.LANES * sim.L * sim.MAX_ITER);

    printf("Average car density: %G\n", density);
    printf("Average flow: %G\n", flow);
    printf("Average lane changes: %G\n", lane_changes);
    printf("Lane changes per car: %f", lane_changes / density);
    printf("Average ping pong lane changes: %G\n", ping_pong_changes);
    
    free(cars);
    free(grid);
    free(new_grid);
    return 0;
}

void initialise(car *cars, int *grid, int *new_grid, parameters sim, int seed) {
    int x, y;
    
    #pragma omp parallel for
    for (int i = 0; i < sim.LANES; i++) {
        for (int j = 0; j < sim.L; j++) {
            grid[sim.L * i + j] = -1;
            new_grid[sim.L * i + j] = -1;
        }
    }

    unsigned short xsubi[3] = {
        seed + 1,
        seed + 2,
        seed + 3
    };

    for (int i = 0; i < sim.N; i++) {
        x = nrand48(xsubi) % sim.L;
        y = nrand48(xsubi) % sim.LANES;
        while (grid[sim.L * y + x] != -1) {
            x = nrand48(xsubi) % sim.L;
            y = nrand48(xsubi) % sim.LANES;
        }
        grid[sim.L * y + x] = i;
        cars[i].x = x;
        cars[i].y = y;
        cars[i].v = 0;
        cars[i].v_d = sim.v_max;
        cars[i].lane_change_now = 0;
    }
}

void update(car *cars, int *grid, int *new_grid, parameters sim, int seed, int iter) {
    #pragma omp parallel for
    for (int i = 0; i < sim.LANES; i++) {
        for (int j = 0; j < sim.L; j++) {
            new_grid[sim.L * i + j] = -1;
        }
    }

    #pragma omp parallel
    {
        unsigned short xsubi[3] = {
            seed + omp_get_thread_num() * iter,
            seed + omp_get_thread_num() * iter + 2,
            seed + omp_get_thread_num() * iter + 3
        };

        #pragma omp for
        for (int i = 0; i < sim.N; i++) {
            cars[i].lane_change_prev = cars[i].lane_change_now;
            if (
                gapAhead(cars, i, grid, sim) < aheadThisLane(cars, i, sim) &&
                gapAheadOther(cars, i, grid, sim) > aheadOtherLane(cars, i, sim) &&
                gapBehindOther(cars, i, grid, sim) > behindOtherLane(cars, i, sim) &&
                erand48(xsubi) < sim.p_change
            ) {
                cars[i].y = (cars[i].y + 1) % sim.LANES;
                cars[i].lane_change_now = 1;
            }
            else {
                cars[i].lane_change_now = 0;
            }
            new_grid[sim.L * cars[i].y + cars[i].x] = i;
        }

        int gap;
        #pragma omp for private(gap)
        for (int i = 0; i < sim.N; i++) {
            gap = gapAhead(cars, i, new_grid, sim);
            if (cars[i].v < cars[i].v_d) {
                cars[i].v++;
            }
            if (cars[i].v > gap) {
                cars[i].v = gap;
            }
            if (cars[i].v > 0 && erand48(xsubi) < sim.p_break) {
                cars[i].v--;
            }
        }
        #pragma omp for
        for (int i = 0; i < sim.LANES; i++) {
            for (int j = 0; j < sim.L; j++) {
                grid[sim.L * i + j] = -1;
            }
        }
        #pragma omp for
        for (int i = 0; i < sim.N; i++) {
            cars[i].x = (cars[i].x + cars[i].v) % sim.L;
            grid[sim.L * cars[i].y + cars[i].x] = i;
        }
    }
}

car * preceedingCar(car *cars, int car_index, int *grid, parameters sim) {
    int lane = cars[car_index].y;
    int x = cars[car_index].x + 1;
    while (grid[sim.L * lane + x] != -1) {
        x = (x + 1) % sim.L;
    }
    return &cars[grid[sim.L * lane + x]];
}

int aheadThisLane(car *cars, int car_index, parameters sim) {
    return cars[car_index].v + 1;
}

int aheadOtherLane(car *cars, int car_index, parameters sim) {
    return cars[car_index].v + 1;
}

int behindOtherLane(car *cars, int car_index, parameters sim) {
    return sim.v_max;
}

int gapAhead(car *cars, int car_index, int *grid, parameters sim) {
    int lane = cars[car_index].y;
    int x = (cars[car_index].x + 1) % sim.L;
    int gap = 0;

    while (grid[sim.L * lane + x] == -1) {
        x = (x + 1) % sim.L;
        gap++;
    }
    return gap;
}

int gapAheadOther(car *cars, int car_index, int *grid, parameters sim) {
    int lane = (cars[car_index].y + 1) % sim.LANES;
    int x = cars[car_index].x;
    int gap = 0;

    if (grid[sim.L * lane + x] != -1) {
        return -1;
    }

    x = (x + 1) % sim.L;
    while (grid[sim.L * lane + x] == -1) {
        x = (x + 1) % sim.L;
        gap++;
    }
    return gap;
}

int gapBehindOther(car *cars, int car_index, int *grid, parameters sim) {
    int lane = (cars[car_index].y + 1) % sim.LANES;
    int x = cars[car_index].x;
    int gap = 0;

    if (grid[sim.L * lane + x] != -1) {
        return -1;
    }

    x = (x - 1) % sim.L;
    while (grid[sim.L * lane + x] == -1) {
        x = mod((x - 1), sim.L);
        gap++;
    }
    return gap;
}

void printGrid(car *cars, int *grid, int length, parameters sim) {
    printf("\033[2A\033[K");
    for (int i = 0; i < sim.LANES; i++) {
        for (int j = 0; j < length; j++) {
            if (grid[sim.L * i + j] != -1) {
                printf("\033[1;32m%d\033[0m", cars[grid[sim.L * i + j]].lane_change_now);
            }
            else {
                printf("\033[2;35m\u2588\033[0m");
            }
        }
        printf("\n\033[K");
    }
}

int mod(int a, int n) {
    return a - n * (int) floor((double) a / (double) n);
}
