#include "solver.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>

/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/

std::vector<bool> create_danger_arr(unsigned int k,
                                    std::vector<unsigned int> partial_soln);

std::vector<std::vector<bool> >
  create_update_arr(unsigned int k,
                        std::vector<unsigned int> partial_soln);

// void recurse_soln(unsigned int start, unsigned int end, unsigned int n,
//                   std::vector<std::vector<unsigned int> > &output_arr,
//                   std::vector<unsigned int> partial_soln,
//                   void (*iter_fn)(std::vector<std::vector<unsigned int> >&),
//                   void (*term_fn)(std::vector<unsigned int>,
//                                   std::vector<std::vector<unsigned int> >&));
void recurse_soln(unsigned int start, unsigned int end, unsigned int n,
                  std::vector<std::vector<unsigned int> > &output_arr,
                  std::vector<unsigned int> partial_soln,
                  std::vector<std::vector<unsigned int> > & queue,
                  int proc_id);

void seq_recursion(
  unsigned int n,
  std::vector<std::vector<unsigned int> >& all_solns,
  std::vector<unsigned int> partialSol,
  unsigned int curN
);

bool check(std::vector<unsigned int>& partialSol);

/*************************** solver.h functions ************************/

// // Serial n-queen solver
// void seq_solver(unsigned int n,
//                 std::vector<std::vector<unsigned int> >& all_solns) {

//   // Initialize
//   std::vector<unsigned int> partial_soln(n, -1);

//   // Create the needed lambda functions
//   auto do_nothing = [](std::vector<std::vector<auto> > &v_arr) {};
//   auto append_result = [](std::vector<auto> v,
//                           std::vector<std::vector<auto> > &v_arr) {
//                          v_arr.push_back(v);
//                        };

//   // Main call
//   recurse_soln(0, n - 1, n, all_solns, partial_soln, do_nothing, append_result);
// }
void seq_solver(
  unsigned int n,
  std::vector<std::vector<unsigned int> >& all_solns
) {
  if (n < 1) {
    return;
  }

  std::vector<unsigned int> v;

  return seq_recursion(n, all_solns, v, 0);
}

// Controls the master processor
void nqueen_master(	unsigned int n,
                    unsigned int k,
                    std::vector<std::vector<unsigned int> >& all_solns) {




	// TODO: Implement this function

	/* Following is a general high level layout that you can follow
	 (you are not obligated to design your solution in this manner.
	  This is provided just for your ease). */

    int flag = 0;
    MPI_STATUS status;
    std::vector<std::vector<unsigned int> > queue;

    // initialize parital sol
    std::vector<unsigned int> partial_soln(-1, n);

    recurse_soln(0, k, n, all_solns, partial_soln, queue, 0);

    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    std::vector<bool> finalCheck(false, num_procs);
    finalCheck[0] = true;

    bool shouldContinue = true;
    while (shouldContinue) {
      flag = 0;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
      if (!flag) {
        // Tag 444: Worker sending soln
        if (status.MPI_TAG == 444) {
          int count = status.count;
          int tmpSz = all_solns.size();
          all_solns.resize(tmpSz + 1);
          all_solns[tmpSz].resize(count);
          MPI_Receive(&all_solns[tmpSz][0], count, MPI_UNSIGNED,
                      status.MPI_SOURCE, 444, MPI_COMM_WORLD,
                      MPI_STATUS_IGNORE);
        } else if (status.MPI_TAG == 999) {
          // Tag 999: Worker sending request for more work
          if (queue.size() == 0) {
            // send kill signal (tag == 4)
            int kill = -1;
            MPI_Send(&kill, 1, MPI_INT, status.MPI_SOURCE, 4, MPI_COMM_WORLD);
            finalCheck[status.MPI_SOURCE] = true;
            shouldContinue = false;
            for (int i = 0; i < finalCheck.size(); i++) {
              if (finalCheck[i] == false) {
                shouldContinue = true;
                break;
              }
            }
          } else {
            // continue
            std::vector<unsigned int> tmp = queue.back();
            queue.pop_back();
            MPI_Send(&tmp, tmp.size(), MPI_UNSIGNED, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
          }
        } else {
          // something is horribly wrong
          // some weird message is being accepted
          std::cout << "Horrible mistake. ABORT\n";
          std::cout << "Error code is just the tag we received hehehe";
          MPI_ABORT(MPI_COMM_WORLD, status.MPI_TAG);
        }



      }
    }
}

// Controls the worker processor
void nqueen_worker(	unsigned int n,
                    unsigned int k) {



	// TODO: Implement this function

	// Following is a general high level layout that you can follow (you are not
  // obligated to design your solution in this manner. This is provided just
  // for your ease).

	/*******************************************************************
	 *
	 * while() {
	 *
	 * 		wait for a message from master
	 *
	 * 		if (message is a partial job) {
	 *				- finish the partial solution
	 *				- send all found solutions to master
	 * 		}
	 *
	 * 		if (message is a kill signal) {
	 *
	 * 				quit
	 *
	 * 		}
	 *	}
	 */

  std::vector<unsigned int> partial_soln;

  // Think of actual end condition later
  while (true) {
    // Not even trying to send data
    // Just send tag
    // Make tag extra(999) so that Master knows to send this boy some work.
    int tempInt = -1;
    MPI_Send(&tempInt, 1, MPI_INT, 0, 999, MPI_COMM_WORLD);

    // first receive size of soln;
    // we have an integer space playing around so let's use that
    MPI_STATUS status;
    MPI_Probe(0, 1, MPI_COMM_WORLD, &status);

    // check for kill signal here
    if (status.MPI_TAG == 4) {
      // welcome to kill signal
      MPI_Receive(&tempInt, 1, MPI_INT, 0, 4, MPI_COMM_WORLD);
      return;
    }

    tempInt = status.count;

    partial_soln.resize(tempInt);
    MPI_Receive(&partial_soln[0], tempInt, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


    // Now we call recurse_soln
    // notice we don't need output_arr NOR do we need queue for worker.
    // start = k, end = n, n = n
    // output_arr = placeholder, partial_soln = partial_soln,
    // queue = placeholder, proc_id = (anything BUT 0);
    std::vector<std::vector<unsigned int> > placeholder;
    recurse_soln(k, n - 1, n, placeholder, partial_soln, placeholder, 1);
  }



}



/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/

/*
Name:
  print_vector

Description:
  simple function used to print a vector
*/
void print_vector(std::vector<auto> v) {
  for (int i = 0; i < v.size(); i++) {
    std::cout << v[i] << " ";
  }
  std::cout << std::endl;
}


/*
Name:
  recurse_soln

Description:
Naive implementation of the above (more flops, less overhead)
*/

void recurse_soln(unsigned int start, unsigned int end, unsigned int n,
                  std::vector<std::vector<unsigned int> > &output_arr,
                  std::vector<unsigned int> partial_soln,
                  std::vector<std::vector<unsigned int> > & queue,
                  int proc_id)
{
                  // void (*iter_fn)(std::vector<std::vector<unsigned int> >&),
                  // void (*term_fn)(std::vector<unsigned int>,
                  //                 std::vector<std::vector<unsigned int> >&)) {
  // Immediately call the iteration function
  // iter_fn(output_arr);
  int flag = 0;
  MPI_STATUS status;
  if (proc_id == 0) {
    flag = 0;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    if (flag && queue.size() != 0) {
      // Tag 1:Worker Sending size of soln [Deprecated] (before being used)
      // Tag 444: Worker sending soln
      if (status.MPI_TAG == 444) {
        int count = status.count;
        int tmpSz = output_arr.size();
        output_arr.resize(tmpSz + 1);
        output_arr[tmpSz].resize(count);
        MPI_Receive(&output_arr[tmpSz][0], count, MPI_UNSIGNED,
                    status.MPI_SOURCE, 444, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
      } else if (status.MPI_TAG == 999) {
        // Tag 999: Worker sending request for more work
        std::vector<unsigned int> tmp = queue.pop_back();
        MPI_Send(&tmp, tmp.size(), MPI_UNSIGNED, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
      } else {
        // something is horribly wrong
        // some weird message is being accepted
        std::cout << "Horrible mistake. ABORT\n";
        std::cout << "Error code is just the tag we received hehehe";
        MPI_ABORT(MPI_COMM_WORLD, status.MPI_TAG);
      }
    }
  }
  // Iter not needed for worker


  // Find all positions to place a queen
  std::vector<unsigned int> q_posns;
  for (int i = 0; i < n; i++) {
    bool safe = true;
    for (int j = 0; j < start; j++) {
      int placed_q = (short) partial_soln[j];
      if (placed_q == i)
        safe = false;
      if ((placed_q - (start - j)) == i)
        safe = false;
      if ((placed_q + (start - j)) == i)
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
    for (int i = 0; i < q_posns.size(); i++) {
      partial_soln[end] = q_posns[i];
      // (*term_fn)(partial_soln, output_arr);
      if (proc_id == 0) {
        queue.push_back(partial_soln);
      } else {
        // we plan to stop this boy here
        int solCount = partial_soln.size();

        // send everythign
        for (int i = 0; i < solCount; i++) {
          MPI_Send(&partial_soln[0], n, MPI_UNSIGNED, 0, 444, MPI_COMM_WORLD);
        }
      }
    }
  }

  // Else, we recurse on the found positions
  else {
    for (int i = 0; i < q_posns.size(); i++) {

      // Assign a copy of the vars
      int q_posn = q_posns[i];
      std::vector<unsigned int> next_partial_soln = partial_soln;

      // Update the partial soln with the incumbent
      next_partial_soln[start] = q_posn;

      // Recurse
      recurse_soln(start + 1, end, n, output_arr, next_partial_soln,
                   queue, proc_id);
    }
  }
  return;
}

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
    if (partialSol[i] == cur || std::abs(size - 1 - i) == std::abs((int)partialSol[i] - (int)cur)) {
      // std::cout << "i:" << i << "|| size:" << size << "|| partialSol[i]" << partialSol[i] << "|| cur:" << cur << std::endl;
      return false;
    }
  }
  return true;
}