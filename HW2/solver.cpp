#include "solver.h"
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>

/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/
void print_vector(std::vector<unsigned int> v);
void recurse_soln(unsigned int start, unsigned int end, unsigned int n,
                  std::vector<std::vector<unsigned int> > &output_arr,
                  std::vector<unsigned int> &partial_soln,
                  std::queue<std::vector<unsigned int> > &queue,
                  int rank);
void seq_recursion(unsigned int n,
                   std::vector<std::vector<unsigned int> >& all_solns,
                   std::vector<unsigned int> partialSol,
                   unsigned int curN);
bool check(std::vector<unsigned int>& partialSol);
/*************************** solver.h functions ************************/

// Serial n-queen solver
void seq_solver(unsigned int n,
                std::vector<std::vector<unsigned int> >& all_solns) {

  // Simple recusive algorithm
  if (n < 1) return;
  std::vector<unsigned int> v;
  return seq_recursion(n, all_solns, v, 0);
}

// Controls the master processor
void nqueen_master(unsigned int n,
                   unsigned int k,
                   std::vector<std::vector<unsigned int> >& all_solns) {

  // --- Main body of the function ---

  // Initialize
  MPI_Status status;
  std::queue<std::vector<unsigned int> > queue;
  std::vector<unsigned int> partial_soln(n, 0);

  // Call the main driver based on use cases
  if (k > 0) {
    std::cout << "Master calls recurse_soln!" << std::endl;
    recurse_soln(0, k - 1, n, all_solns, partial_soln, queue, 0);
  } else {
    std::cout << "Send all the work to one processor!" << std::endl;
    // Acknowledge the request
    int tag_nmbr;
    MPI_Status status;
    MPI_Recv(&tag_nmbr, 1, MPI_INT, MPI_ANY_SOURCE, 999, MPI_COMM_WORLD,
             &status);
    // Now send an empty solution to the lucky processor
    std::cout << "Master sends (only once): ";
    print_vector(partial_soln);
    MPI_Send(&partial_soln[0], n, MPI_UNSIGNED, status.MPI_SOURCE, 1,
             MPI_COMM_WORLD);
  }

  // Begin termination phase
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  std::vector<bool> finalCheck(num_procs, false);
  finalCheck[0] = true;
  bool shouldContinue = true;
  std::vector<unsigned int> complete_soln(n, 0);
  std::cout << std::endl << "Starting termination..." << std::endl;
  std::cout << "Current solutions: " << all_solns.size() <<  std::endl;
  while (shouldContinue) {

    // Check current conditions
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    std::cout << "Tag " << status.MPI_TAG <<
      " received from worker " << status.MPI_SOURCE << std::endl;

    // Tag 444: Worker sending soln
    if (status.MPI_TAG == 444) {
      // Receive it and append it
      MPI_Recv(&complete_soln[0], n, MPI_UNSIGNED, status.MPI_SOURCE,
               444, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::cout << "Master receives: ";
      print_vector(complete_soln);
      all_solns.push_back(complete_soln);
    }

    // Tag 999: Worker sending request for more work
    else if (status.MPI_TAG == 999) {
      // No more jobs available
      if (queue.size() == 0) {
        std::cout << "I'm gonna kill worker " << status.MPI_SOURCE <<
          " now..." << std::endl;
        // send kill signal (tag = 4)
        int kill = -1;
        MPI_Send(&kill, 1, MPI_INT, status.MPI_SOURCE, 4, MPI_COMM_WORLD);
        finalCheck[status.MPI_SOURCE] = true;
        shouldContinue = false;
        for (unsigned int i = 0; i < finalCheck.size(); i++) {
          if (finalCheck[i] == false) {
            shouldContinue = true;
            break;
          }
        }
      }
      // There are still some jobs remaining
      else {
        // Acknowledge the request
        int tag_nmbr;
        MPI_Recv(&tag_nmbr, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Now send out a task
        std::vector<unsigned int> incumbent = queue.front();
        queue.pop();
        MPI_Send(&incumbent[0], n, MPI_UNSIGNED, status.MPI_SOURCE, 1,
                 MPI_COMM_WORLD);
      }
    }

    // something is horribly wrong
    // some weird message is being accepted
    else {
      std::cout << "Horrible mistake. ABORT\n";
      std::cout << "Error code is just the tag we received hehehe";
      MPI_Abort(MPI_COMM_WORLD, status.MPI_TAG);
    }
  }
  // End master's job
  std::cout << std::endl << "Total solutions = " << all_solns.size() <<
    std::endl;

}

// Controls the worker processor
void nqueen_worker(unsigned int n,
                   unsigned int k) {

  // --- Main body of the function ---

  // Initialize
  std::vector<unsigned int> partial_soln(n, 0);
  std::vector<std::vector<unsigned int> > empty_vvec;
  std::queue<std::vector<unsigned int> > empty_queue;
  int proc_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

  // Think of actual end condition later
  while (true) {
    // Not even trying to send data
    // Just send tag
    // Make tag = 999 so that master knows to send this boy some work.
    int tempInt = 999;
    std::cout << "Worker is waiting..." << std::endl;
    MPI_Send(&tempInt, 1, MPI_INT, 0, 999, MPI_COMM_WORLD);

    // first receive size of soln;
    // we have an integer space playing around so let's use that
    MPI_Status status;
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    // check for kill signal here
    if (status.MPI_TAG == 4) {
      // welcome to kill signal
      MPI_Recv(&tempInt, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::cout << "Worker number " << proc_id
                << " is now dying. What a cruel world..." << std::endl;
      return;
    }

    // Master has sent you work
    else if (status.MPI_TAG == 1) {
      // Receive it
      MPI_Recv(&partial_soln[0], n, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      std::cout << "Worker receives: ";
      print_vector(partial_soln);
    }

    // something is horribly wrong
    // some weird message is being accepted
    else {
      std::cout << "Horrible mistake. ABORT\n";
      std::cout << "Error code is just the tag we received hehehe";
      MPI_Abort(MPI_COMM_WORLD, status.MPI_TAG);
    }

    // Now we call recurse_soln
    // notice we don't need output_arr NOR do we need queue for worker.
    // start = k, end = n - 1, n = n
    // output_arr = empty_vvec, partial_soln = partial_soln,
    // queue = empty_queue;
    std::vector<std::vector<unsigned int> > placeholder;
    std::cout << "Worker starts looking for solutions..." << std::endl;
    recurse_soln(k, n - 1, n, empty_vvec, partial_soln, empty_queue, 1);
    std::cout << "Worker stops looking for solutions..." << std::endl;
  }
  // End worker's job

}



/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/

/*
Name:
  print_vector

Description:
  Simple vector printing function for vector<unsigned int> objects
*/
void print_vector(std::vector<unsigned int> v) {
  for (unsigned int i = 0; i < v.size(); i++) {
    std::cout << v[i] << " ";
  }
  std::cout << std::endl;
}

/*
Name:
  recurse_soln

Description:
  Main work horse of the master and worker processors. Finds partial or full
  solutions and communicates these in a (master <-> worker) network, where there
  is only one master and many workers.
*/

void recurse_soln(unsigned int start, unsigned int end, unsigned int n,
                  std::vector<std::vector<unsigned int> > &output_arr,
                  std::vector<unsigned int> &partial_soln,
                  std::queue<std::vector<unsigned int> > &queue,
                  int rank) {

  // Find all safe positions to place a queen
  std::vector<unsigned int> q_posns;
  int start_idx = (short) start;
  int n_idx = (short) n;
  for (int i = 0; i < n_idx; i++) {
    bool safe = true;
    for (int j = 0; j < start_idx; j++) {
      int placed_q = (short) partial_soln[j];
      if (placed_q == i)
        safe = false;
      if ((placed_q - (start_idx - j)) == i)
        safe = false;
      if ((placed_q + (start_idx - j)) == i)
        safe = false;
    }
    if (safe)
      q_posns.push_back(i);
  }

  // If no positions are found, do nothing
  if (q_posns.size() == 0) {
    // Nothing to see here; move along...
  }

  // If we are at the last needed column, call our termination function on
  // the set of completed solutions
  else if (start == end) {
    int q_size = (short) q_posns.size();
    for (int i = 0; i < q_size; i++) {
      partial_soln[end] = q_posns[i];
      if (rank == 0 && end < n - 1) {
        std::cout << "Master creates: ";
        print_vector(partial_soln);
        queue.push(partial_soln);
      } else if (rank == 0 && end == n-1) {
        output_arr.push_back(partial_soln);
      } else {
        std::cout << "Worker is sending a solution: ";
        print_vector(partial_soln);
        MPI_Send(&partial_soln[0], n, MPI_UNSIGNED, 0, 444, MPI_COMM_WORLD);
      }
    }
  }

  // Else, we recurse on the found positions
  else {
    int q_size = (short) q_posns.size();
    for (int i = 0; i < q_size; i++) {
      // Update the partial soln with the incumbent
      partial_soln[start] = q_posns[i];
      // Recurse
      recurse_soln(start + 1, end, n, output_arr, partial_soln, queue, rank);
    }
  }

  // If you are Master, do some checks EVERY call of compute_soln
  if (rank == 0) {
    MPI_Status status;
    int flag = 0;
    // Do this in async
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    // Do something only if Iprobe receives a message; otherwise do nothing
    if (flag) {
      std::cout << "Master sees: " << "queue.size() = " << queue.size() <<
        ", tag = " << status.MPI_TAG << std::endl;

      // Tag 444: Worker sending soln
      if (status.MPI_TAG == 444) {
        // Receive once and continue the computation
        std::vector<unsigned int> complete_soln(n, -1);
        MPI_Recv(&complete_soln[0], n, MPI_UNSIGNED,
                 status.MPI_SOURCE, 444, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        std::cout << "Master receives: ";
        print_vector(complete_soln);
        output_arr.push_back(complete_soln);
      }

      // Tag 999: Worker sending request for more work
      else if (status.MPI_TAG == 999) {
        if (queue.size() != 0) {
          // Acknowledge the request
          int tag_nmbr;
          MPI_Recv(&tag_nmbr, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          // Now send out some tasks
          std::vector<unsigned int> incumbent = queue.front();
          queue.pop();
          std::cout << "Master sends: ";
          print_vector(incumbent);
          MPI_Send(&incumbent[0], n, MPI_UNSIGNED, status.MPI_SOURCE, 1,
                   MPI_COMM_WORLD);
        }
      }

      // something is horribly wrong
      // some weird message is being accepted
      else {
        std::cout << "Horrible mistake. ABORT\n";
        std::cout << "Error code is just the tag we received hehehe";
        MPI_Abort(MPI_COMM_WORLD, status.MPI_TAG);
      }
    }
  }

  // End your call and go up the stack
  return;
}

// For seq_solver ONLY
void seq_recursion(unsigned int n,
          std::vector<std::vector<unsigned int> >& all_solns,
          std::vector<unsigned int> partialSol,
          unsigned int curN)
{
  if (n == 0) {
    // no way this gets activated but just in case.
    std::cout << "Scream and shout! The End Is Coming!";
    return;
  }
  if (curN < n) {
    for (int next = 0; next < (int)n; next++) {
      partialSol.push_back(next);
      if (check(partialSol)) {
        seq_recursion(n, all_solns, partialSol, curN + 1);
      }
      partialSol.pop_back();
    }
  }
  if (curN == n) {
    // [END condition]
    all_solns.push_back(partialSol);
  }
}

bool check(std::vector<unsigned int>& partialSol) {
  int size = partialSol.size();
  unsigned int cur = partialSol[size - 1];
  for (int i = 0; i < size - 1; i++) {
    if (partialSol[i] == cur ||
        std::abs(size - 1 - i) == std::abs((int)partialSol[i] - (int)cur)) {
      return false;
    }
  }
  return true;
}
