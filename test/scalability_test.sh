#!/bin/bash

if [ $# -eq 0 ]; then
	REPS=1
else
	REPS=$1
fi

REP_FILE="./cache.dat"
PREV_REPS=$(cat "$REP_FILE")

cd ../src/
FILE=./performance_test
if [[ $REPS != $PREV_REPS || ! -f $FILE ]]; then
	echo "Building Test Executable"
	if [[ -f $FILE.o ]]; then
		rm $FILE.o
	fi
	make -s performance_test REPS=$REPS
fi

for i in 1 2 4
do
	echo "Test with $i processors"
	mpiexec -np $i ./performance_test
done

cd ../test
echo "$REPS" > "$REP_FILE"
python3 plot_results.py
