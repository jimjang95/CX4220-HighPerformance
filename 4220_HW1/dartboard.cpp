#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>

// using namespace std;

int dboard(int N);

int main(int argc, char *argv[])
{
    if(argc != 3) {
        //TODO: make this useful
        return 0;
    }

    int N = atoi(argv[1]);
    int R = atoi(argv[2]);

    if(R > 100 || N < 0) {
        //TODO: change N check && make useful
        return 0;
    }

    int tmp1;
    double tmp2;
    double avg = 0;
    for(int i = 0; i < R; i++) {
        tmp1 = dboard(N);
        tmp2 = tmp1 / (double)R;
        avg += tmp2;
        // std::cout << "dboard(N), /100, curAvg || " << tmp1 << ' ' << tmp2 << ' ' << avg << std::endl;
        // avg += dboard(N) / (double)R;
    }

    std::ofstream fout;
    fout.open("output.txt");
    fout << "N, R, average: " << N << ", " << R  << ", " << avg << '\n';
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

    //TODO: later change the seed to taskID(rank)
    //[BUG] TODO: this somehow gives the exact same seed everytime. fix it.
    std::srand(std::time(NULL));

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