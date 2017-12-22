=======================================================
How to perform full test of i-PI QuantumEspresso driver
=======================================================

Cases when G-Vectors need to be recalculated are not easily achievable.
However, we can decrease "gvec_omega_tol" parameter to make it happen.

1. Change $QE_ROOT/PW/src/run_driver.f90, line 37:
`REAL*8, PARAMETER :: gvec_omega_tol= 1.0D-1`
to 
`REAL*8, PARAMETER :: gvec_omega_tol= 1.0D-3`

2. Recompile Q-E

3. Go to `$IPI_ROOT/examples/regtest/qespresso`:
`cd $IPI_ROOT/examples/regtest/qespresso`

4. Create references for regtest. Don't try to modify this command, as
regtest is broken and it will probably fail.
`../../../tools/py/regtest.py --tests-folder . --folder-run run --create-reference`

5. Introduce changes to Q-E driver (recompile the code to new version,
link pw.x to different version or do any other relevant thing)

6. Run the regtest (again, do not modify this command):
`../../../tools/py/regtest.py --tests-folder . --folder-run run`

7. There might be numerical discrepancies. To examine differences, compare files in
`relevant_dir/regtest_ref` with files in `run/relevant_dir/io`
