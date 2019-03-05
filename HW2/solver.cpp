#include "solver.h"
#include <cmath>


/*************************** DECLARE YOUR HELPER FUNCTIONS HERE ************************/
void seq_recursion(
	unsigned int n,
	std::vector<std::vector<unsigned int> >& all_solns,
	std::vector<unsigned int> partialSol,
	unsigned int curN
);

bool check(std::vector<unsigned int>& partialSol);


/*************************** solver.h functions ************************/


void seq_solver(
	unsigned int n,
	std::vector<std::vector<unsigned int> >& all_solns
) {

	// TODO: Implement this function
	if (n < 1) {
		return;
	}

	std::vector<unsigned int> v;
	// for (unsigned int i = 0; i < n; i++) {
	// 	v.push_back(i);
	// }
	// for (std::vector<unsigned int>::const_iterator i = v.begin(); i != v.end(); ++i) {

	// 	std::cout << *i << ' ';

	// }
	return seq_recursion(n, all_solns, v, 0);
}






void nqueen_master(	unsigned int n,
					unsigned int k,
					std::vector<std::vector<unsigned int> >& all_solns) {




	// TODO: Implement this function

	/* Following is a general high level layout that you can follow
	 (you are not obligated to design your solution in this manner.
	  This is provided just for your ease). */


	/******************* STEP 1: Send one partial solution to each worker ********************/
	/*
	 * for (all workers) {
	 * 		- create a partial solution.
	 * 		- send that partial solution to a worker
	 * }
	 */


	/******************* STEP 2: Send partial solutions to workers as they respond ********************/
	/*
	 * while() {
	 * 		- receive completed work from a worker processor.
	 * 		- create a partial solution
	 * 		- send that partial solution to the worker that responded
	 * 		- Break when no more partial solutions exist and all workers have responded with jobs handed to them
	 * }
	 */

	/********************** STEP 3: Terminate **************************
	 *
	 * Send a termination/kill signal to all workers.
	 *
	 */





}

void nqueen_worker(	unsigned int n,
					unsigned int k) {



	// TODO: Implement this function

	// Following is a general high level layout that you can follow (you are not obligated to design your solution in this manner. This is provided just for your ease).

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


}



/*************************** DEFINE YOUR HELPER FUNCTIONS HERE ************************/

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




