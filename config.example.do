/* Copy to config.local.do; never commit private keys. */
global H_CENSUS_KEY "INSERT_YOUR_CENSUS_KEY"
/* mirror = documented third-party BLS LAUS distributor, automatic.
   official = original BLS flat files, which may return HTTP 403. */
global H_LAUS_SOURCE_MODE "mirror"
/* Definitive four-outcome placebo refits. 199 draws by default. */
global H_JOINT_PLACEBO_REPS 199
global H_JOINT_PLACEBO_SEED 20261114
global H_JOINT_PLACEBO_ITER 30000
global H_JOINT_PLACEBO_TOL 0.00005
global H_JOINT_PLACEBO_STRICT 1
global H_RUN_CENSUS 1
global H_CENSUS_FILE ""
/* Optional: full path to a python3 that has NumPy, pandas and Matplotlib.
   Leave empty to search python3, /usr/local/bin/python3, Homebrew and
   python.org locations in that order. Stata's shell PATH does not include
   /usr/local/bin, so a bare `python3` there is Apple's dependency-free one. */
global H_PYTHON ""
