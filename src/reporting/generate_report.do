/* Readable report and full graph atlas, gated on current-run stage statuses.
   STALE FILE IGNORED: the Python renderer excludes every non-PASSED family.
   Numerical tables are read only; no estimator is rerun here.

   Revision 14 Nov 2026: Stata's `shell` runs with PATH=/usr/bin:/bin:/usr/sbin:
   /sbin, so on macOS a bare `python3` resolves to Apple's CommandLineTools
   interpreter, which has no NumPy, even when the Python the user installs
   packages into is on their own PATH. Installing the requirements from a
   terminal therefore does NOT fix this stage. Try the usual interpreter
   locations in turn and keep the first that builds the report; set
   global H_PYTHON in config.local.do to force a specific one. */
version 16.0
local root "${H_ROOT}"
if "`root'"=="" local root "`c(pwd)'"
local okfile "`root'/outputs/audit/REPORT_BUILD_OK.txt"
local cands `""${H_PYTHON}" "python3" "/usr/local/bin/python3""'
local cands `"`cands' "/opt/homebrew/bin/python3""'
local cands `"`cands' "/Library/Frameworks/Python.framework/Versions/Current/bin/python3""'
local cands `"`cands' "/usr/bin/python3""'
local used ""
local tried ""
foreach p of local cands {
    if "`p'"=="" continue
    if "`used'"!="" continue
    capture erase "`okfile'"
    capture shell "`p'" "`root'/src/reporting/make_report.py" --root "`root'"
    capture confirm file "`okfile'"
    if !_rc local used "`p'"
    else local tried `"`tried' `p'"'
}
if "`used'"=="" {
    di as error "Reporting failed with every Python tried:`tried'"
    di as error "Install the dependencies in ONE of them:"
    di as error "  <python3> -m pip install -r src/reporting/requirements.txt"
    di as error `"then put its full path in config.local.do: global H_PYTHON "/full/path/to/python3""'
    di as error "Rerun run_report_only.do; no estimate is lost."
    exit 601
}
if `"`tried'"'!="" di as txt "Python tried and rejected (missing NumPy/pandas/Matplotlib):`tried'"
di as result "Report built with `used'"
di as result "Readable results: outputs/reports/HELENE_COMPLETE_RESULTS.txt"
di as result "Graph atlas: outputs/reports/HELENE_GRAPH_ATLAS.pdf"
di as result "Graph gallery: outputs/reports/graphs/index.html"
