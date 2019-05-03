#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
    int LANE;           // Number of lanes
    int L;              // Length of the road
    int N;              // Number of cars
    int MAX_ITER;       // Maximum iterations
    int v_max;          // Maximum desired speed
    double p_change;    // Lane change probability
    double p_break;     // Random breaking probability
} simulation_params;

void initialize(int *road, int *v, simulation_params *sim_params, int seed);
void update(int *road, int *interm_road, int *v, int *interm_v, simulation_params *sim_params, unsigned short *xsubi);
int aheadParam(int v, simulation_params *sim_params);
int aheadOtherParam(int v, simulation_params *sim_params);
int behindOtherParam(int v, simulation_params *sim_params);
int gapAhead(int lane, int position, int *road, simulation_params *sim_params);
int gapAheadOther(int lane, int position, int *road, simulation_params *sim_params);
int gapBehindOther(int lane, int position, int *road, simulation_params *sim_params);

void printRoad(int *road, int LANE, int L);
int mod(int a, int n);

int main(int argc, char *argv[]) {
    if (argc == 1) {
        printf("Usage: ./traffic random_seed\n");
        printf("random_seed: integer\n");
        exit(99);
    }

    // Set simulation parameters
    simulation_params sim_params = {
        2,
        133333,
        sim_params.LANE * sim_params.L * 0.08,
        5000,
        5,
        1,
        0.5
    };

    struct timespec sleeper = {     // Sleep for visualization
        0,
        300000000
    };
    struct timespec remain;

    printf("\n\n");     // Creating space to print out the cars

    int *road = malloc(sim_params.LANE * sim_params.L * sizeof(int));         // Main road grid
    int *interm_road = malloc(sim_params.LANE * sim_params.L * sizeof(int));  // Road grid for intermediate calculations
    int *v = malloc(sim_params.LANE * sim_params.L * sizeof(int));            // Speed grid
    int *interm_v = malloc(sim_params.LANE * sim_params.L * sizeof(int));     // Speed grid for intermediate calculations

    // Statistics to be collected
    double density = (double) sim_params.N / (double) (sim_params.LANE * sim_params.L);
    double flow = 0;
    
    int seed = atoi(argv[1]);                           // Random seed provided at program invoke
    
    if (road == NULL || v == NULL) {
        printf("Could not allocate memory.\n");
        exit(1);
    }
    initialize(road, v, &sim_params, seed);

    // printf("Original Road:\n\n");
    // printRoad(road, LANE, 10);
    // printf("\n");
    
    unsigned short xsubi[3];
    xsubi[0] = seed + 1;
    xsubi[1] = seed + 2;
    xsubi[2] = seed + 3;
    for (int i = 0; i < 1000; i++) {
        for (int lane = 0; lane < sim_params.LANE; lane++) {
            for (int j = 0; j < sim_params.L; j++) {
                interm_road[sim_params.L*lane + j] = 0;
                interm_v[sim_params.L*lane + j] = 0;
            }
        }
        for (int lane = 0; lane < sim_params.LANE; lane++) {
            for (int j = 0; j < sim_params.L; j++) {
                if (road[sim_params.L*lane + j] == 1) {
                    if (
                        gapAhead(lane, j, road, &sim_params) < aheadParam(v[sim_params.L*lane + j], &sim_params) &&
                        gapAheadOther(lane, j, road, &sim_params) > aheadOtherParam(v[sim_params.L*lane + j], &sim_params) &&
                        gapBehindOther(lane, j, road, &sim_params) > behindOtherParam(v[sim_params.L*lane + j], &sim_params) &&
                        erand48(xsubi) < sim_params.p_change
                    ) {
                        interm_road[sim_params.L*mod(lane + 1, sim_params.LANE) + j] = 1;
                        interm_v[sim_params.L*mod(lane + 1, sim_params.LANE) + j] = v[sim_params.L*lane + j];
                    }
                    else {
                        interm_road[sim_params.L*lane + j] = 1;
                        interm_v[sim_params.L*lane + j] = v[sim_params.L*lane + j];
                    }
                }
            }
        }
        for (int lane = 0; lane < sim_params.LANE; lane++) {
            for (int j = 0; j < sim_params.L; j++) {
                road[sim_params.L*lane + j] = 0;
                v[sim_params.L*lane + j] = 0;
            }
        }
        for (int lane = 0; lane < sim_params.LANE; lane++) {
            int gap;
            for (int j = 0; j < sim_params.L; j++) {
                if (interm_road[sim_params.L*lane + j] == 1) {
                    gap = gapAhead(lane, j, interm_road, &sim_params);
                    if (interm_v[sim_params.L*lane + j] != sim_params.v_max) {
                        interm_v[sim_params.L*lane + j]++;
                    }
                    if (interm_v[sim_params.L*lane + j] > gap) {
                        interm_v[sim_params.L*lane + j] = gap;
                    }
                    if (interm_v[sim_params.L*lane + j] != 0 && erand48(xsubi) < sim_params.p_break) {
                        interm_v[sim_params.L*lane + j] -= 1;
                    }
                    road[sim_params.L*lane + mod(j + interm_v[sim_params.L*lane + j], sim_params.L)] = 1;
                    v[sim_params.L*lane + mod(j + interm_v[sim_params.L*lane + j], sim_params.L)] = interm_v[sim_params.L*lane + j];
                }
            }
        }
        printRoad(road, sim_params.LANE, 100);
        nanosleep(&sleeper, &remain);
    }

    for (int i = 0; i < sim_params.MAX_ITER; i++) {
        update(road, interm_road, v, interm_v, &sim_params, xsubi);
        if (i % 5 == 0) {
            for (int lane = 0; lane < sim_params.LANE; lane++) {
                for (int j = 0; j < sim_params.L; j++) {
                    flow += v[sim_params.L*lane + j];
                }
            }
        }
    }
    flow /= (double) (sim_params.L * sim_params.LANE * sim_params.MAX_ITER);
    // int cars = 0;
    // for (int lane = 0; lane < LANE; lane++) {
    //     for (int j = 0; j < L; j++) {
    //         cars += road[L*lane + j];
    //     }
    // }
    // printf("%d\n", cars);
    // printf("Road after %d iterations:\n\n", MAX_ITER);
    // printRoad(road, LANE, 10);
    // printf("\n");

    printf("Average road density: %f\n", density);
    printf("Traffic flow (time and space averaged speed): %f\n", flow);
    
    free(road);
    free(interm_road);
    free(v);
    free(interm_v);
    return 0;
}

void initialize(int *road, int *v, simulation_params *sim_params, int seed) {
    int x, y;
    
    for (int i = 0; i < sim_params->LANE * sim_params->L; i++) {
        road[i] = 0;
        v[i] = 0;
    }

    unsigned short xsubi[3];
    xsubi[0] = seed + 1;
    xsubi[1] = seed + 2;
    xsubi[2] = seed + 3;
    
    x = nrand48(xsubi) % sim_params->LANE;
    y = nrand48(xsubi) % sim_params->L;
    for (int i = 0; i < sim_params->N; i++) {
        road[sim_params->L*x + y] = 1;
        while (road[sim_params->L*x + y] == 1) {
            x = nrand48(xsubi) % sim_params->LANE;
            y = nrand48(xsubi) % sim_params->L;
        }
    }
}

void update(int *road, int *interm_road, int *v, int *interm_v, simulation_params *sim_params, unsigned short *xsubi) {
    for (int lane = 0; lane < sim_params->LANE; lane++) {
        for (int j = 0; j < sim_params->L; j++) {
            interm_road[sim_params->L*lane + j] = 0;
            interm_v[sim_params->L*lane + j] = 0;
        }
    }
    for (int lane = 0; lane < sim_params->LANE; lane++) {
        for (int j = 0; j < sim_params->L; j++) {
            if (road[sim_params->L*lane + j] == 1) {
                if (
                    gapAhead(lane, j, road, sim_params) < aheadParam(v[sim_params->L*lane + j], sim_params) &&
                    gapAheadOther(lane, j, road, sim_params) > aheadOtherParam(v[sim_params->L*lane + j], sim_params) &&
                    gapBehindOther(lane, j, road, sim_params) > behindOtherParam(v[sim_params->L*lane + j], sim_params) &&
                    erand48(xsubi) < sim_params->p_change
                ) {
                    interm_road[sim_params->L*mod(lane + 1, sim_params->LANE) + j] = 1;
                    interm_v[sim_params->L*mod(lane + 1, sim_params->LANE) + j] = v[sim_params->L*lane + j];
                }
                else {
                    interm_road[sim_params->L*lane + j] = 1;
                    interm_v[sim_params->L*lane + j] = v[sim_params->L*lane + j];
                }
            }
        }
    }
    for (int lane = 0; lane < sim_params->LANE; lane++) {
        for (int j = 0; j < sim_params->L; j++) {
            road[sim_params->L*lane + j] = 0;
            v[sim_params->L*lane + j] = 0;
        }
    }
    for (int lane = 0; lane < sim_params->LANE; lane++) {
        int gap;
        for (int j = 0; j < sim_params->L; j++) {
            if (interm_road[sim_params->L*lane + j] == 1) {
                gap = gapAhead(lane, j, interm_road, sim_params);
                if (interm_v[sim_params->L*lane + j] != sim_params->v_max) {
                    interm_v[sim_params->L*lane + j]++;
                }
                if (interm_v[sim_params->L*lane + j] > gap) {
                    interm_v[sim_params->L*lane + j] = gap;
                }
                if (interm_v[sim_params->L*lane + j] != 0 && erand48(xsubi) < sim_params->p_break) {
                    interm_v[sim_params->L*lane + j] -= 1;
                }
                road[sim_params->L*lane + mod(j + interm_v[sim_params->L*lane + j], sim_params->L)] = 1;
                v[sim_params->L*lane + mod(j + interm_v[sim_params->L*lane + j], sim_params->L)] = interm_v[sim_params->L*lane + j];
            }
        }
    }
}

int aheadParam(int v, simulation_params *sim_params) {
    return v + 1;
}

int aheadOtherParam(int v, simulation_params *sim_params) {
    return v + 1;
}

int behindOtherParam(int v, simulation_params *sim_params) {
    return sim_params->v_max;
}

int gapAhead(int lane, int position, int *road, simulation_params *sim_params) {
    int j = mod(position + 1, sim_params->L);
    int count = 0;
    while (road[sim_params->L*lane + j] != 1) {
        j = mod(j + 1, sim_params->L);
        count++;
        if (j == position) {
            break;
        }
    }
    return count;
}

int gapAheadOther(int lane, int position, int *road, simulation_params *sim_params) {
    int j = mod(position + 1, sim_params->L);
    int count = 0;
    lane = mod(lane + 1, sim_params->LANE);

    if (road[sim_params->L*lane + position] == 1) {
        return -1;
    }

    while (road[sim_params->L*lane + j] != 1) {
        j = mod(j + 1, sim_params->L);
        count++;
        if (j == position) {
            break;
        }
    }
    return count;
}

int gapBehindOther(int lane, int position, int *road, simulation_params *sim_params) {
    int j = mod(position - 1, sim_params->L);
    int count = 0;
    lane = mod(lane + 1, sim_params->LANE);

    if (road[sim_params->L*lane + position] == 1) {
        return -1;
    }

    while (road[sim_params->L*lane + j] != 1) {
        j = mod(j + 1, sim_params->L);
        count++;
        if (j == position) {
            break;
        }
    }
    return count;
}

void printRoad(int *road, int LANE, int L) {
    printf("\033[2A\033[K");
    for (int i = 0; i < LANE; i++) {
        for(int j = 0; j < L; j++) {
            if (road[L*i + j] == 1) {
                printf("\033[01;32m\u2588\033[00m");
            }
            else {
                printf("\033[02;35m\u2588\033[00m");
            }
        }
        printf("\n");
    }
}

int mod(int a, int n) {
    return a - n * (int) floor((double) a / (double) n);
}
