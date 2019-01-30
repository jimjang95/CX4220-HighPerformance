EXECS=prog1
MPICC?=mpicxx

all: ${EXECS}

prog1: prog1.cpp
	${MPICC} -std=c++11 prog1.cpp -o prog1

clean:
	rm -f ${EXECS} output.txt