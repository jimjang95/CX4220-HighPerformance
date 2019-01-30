#include <iostream>
#include <mpi.h>
#include <cmath>
#include <string>

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

    int N; // ready to receive from Master
    int R; // ready to receive from Master
    double pi = 0; // to be used by Master

    // Step 1: read / broadcast 'N'
    if (rank == 0) {
        // MASTER
        // init checks
        if (argc != 3) {
            std::cerr << "Not enough arguments\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
            MPI_Finalize();
            return 1;
        }
        N = std::stoi(argv[1]);
        R = std::stoi(argv[2]);

        // init checks
        if (N < 5000000) {
            std::cerr << "N is an invalid value. Input 5000000 < N.\n";
            MPI_Abort(MPI_COMM_WORLD, 2);
            MPI_Finalize();
            return 2;
        }
        if (R > 100 || R < 1) {
            std::cerr << "R is an invalid value. Input 0 < R < 100.\n";
            MPI_Abort(MPI_COMM_WORLD, 2);
            MPI_Finalize();
            return 2;
        }
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&R, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int part;
    int total;
    for (int i = 0; i < R; i++) {
        // Step 2: everyone does dboard
        part = dboard(N / worldSize);

        // Step 3: add together all the dboard values.
        MPI_Reduce(&part, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            pi += (double)N * 2 / total;
        }
    }

    if (rank == 0) {
        pi = pi / R;
        time += MPI_Wtime();
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
    double square = std::sqrt(r / 2.0f);

    double a;
    double theta;

    int M = 0;

    for(int i = 0; i < N; i++) {
        a = (double)std::rand() / RAND_MAX * r;
        a = std::sqrt(a);
        theta = (double)std::rand() / RAND_MAX * 360;
        if(std::abs(a * std::sin(theta)) <= square \
            && std::abs(a * std::cos(theta)) <= square)
        {
            M += 1;
        }
    }
    return M;
}