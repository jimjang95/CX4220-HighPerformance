#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <string>
// #include <fstream>

int dboard(int N);

int main(int argc, char *argv[])
{
    // set up MPI
    MPI_Init(&argc, &argv);
    double time = -1 * MPI_Wtime();

    // get # of processes
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    // this is the rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::srand(rank); // seed here
    // printf("[RAND CHECK] 1");
    // printf("Rank %d/%d: ", rank, worldSize);
    // printf("%d, %d, %d\n", std::rand(), std::rand(), std::rand());

    int N; // ready to receive from Master
    int R; // ready to receive from Master
    double pi = 0; // to be used by Master

    // Step 1: read / broadcast 'N'
    if (rank == 0) {
        // MASTER
        N = atoi(argv[1]);
        // R = atoi(argv[2]);
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&R, 1, MPI_INT, 0, MPI_COMM_WORLD);
    R = atoi(argv[2]);
    // printf("\n");

    int part;
    int total;
    for (int i = 0; i < R; i++) {
        // Step 2: everyone does dboard
        // printf("Rank %d/%d: ", rank, worldSize);
        // printf("N: %d , N/worldSize: %d\n", N, N / worldSize);
        part = dboard(N / worldSize);
        // double part = (double)M / N;
        // printf("Rank %d/%d: Sending in %d shots. ", rank, worldSize, part);
        // printf("part = %f.\n", part);

        // Step 3: add together all the dboard values.
        // MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&part, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            pi += (double)N * 2 / total;
            // printf("pi: %f\n", pi);
        }
        // printf("Rank %d/%d: ", rank, worldSize);
        // printf("total, pi = %d, %f.\n", total, pi);
    }

    if (rank == 0) {
        // printf("pi, R = %f, %d.\n", pi, R);
        pi = pi / R;
        time += MPI_Wtime();
        // printf("Rank %d/%d: ", rank, worldSize);
        // printf("pi = %f.\n", pi);
        std::string output = "";
        output += "N= ";
        output += std::to_string(N);
        output += ", R= ";
        output += std::to_string(R);
        output += ", P= ";
        output += std::to_string(worldSize);
        output += ", PI= ";
        output += std::to_string(pi);
        output += "\nTime = ";
        output += std::to_string(time);
        output += "s\n";
        std::cout << output;
    }

    // ze end
    MPI_Finalize();
    return 0;
}

// int N: number of darts to throw
int dboard(int N) {
    // Basic board size set
    double r = 2.0;
    double square = r / std::sqrt(2.0);

    double a;
    double theta;

    int M = 0;

    for(int i = 0; i < N; i++) {
        a = (double)std::rand() / RAND_MAX * r;
        theta = (double)std::rand() / RAND_MAX * 360;
        if(std::abs(a * std::sin(theta)) <= square \
            && std::abs(a * std::cos(theta)) <= square)
        {
            // std::cout << "a: " << a << '\n';
            // std::cout << "a sin: " << std::abs(a * std::sin(theta)) << '\n';
            // std::cout << "a cos: " << std::abs(a * std::cos(theta)) << '\n';
            // std::cout << "square: " << square << '\n';

            M += 1;
        }
    }
    return M;
}