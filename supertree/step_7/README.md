# Step 7 Decompose the inputs into subproblems of uncontested taxa
Most of the work is actually done in a scratch directory - see [that README.md](../step_7_scratch/README.md)

Once the subproblems are obtained and check-summed, `make` will execute:
    python move-subproblems-if-differing.py step_7_scratch/checksummed-subproblem-ids.txt step_7_scratch/export-sub-temp step_7 step_8 step_9

This requires [peyotl](https://github.com/OpenTreeOfLife/peyotl) code that is currently on the 
"supertree" branch of that repository. The `move-subproblems-if-differing.py` tool:
  1. treats the first argument, `step_7_scratch/checksummed-subproblem-ids.txt`, as the complete list of valid subproblems
  2. treats the second argument as the parent of the dumped, check-summed subproblems.
  3. writes any changed subproblem to `step_7`
  4. if the subproblem is trivially solvable, writes a solution to `step_9`, otherwise it writes a simplification of the
    subproblem to `step_8`.