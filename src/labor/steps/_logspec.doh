/* _logspec.doh -- append one estimated specification to the regression log.

   Included (NOT run) immediately after an estimation command, so it shares the
   caller's local macros and e() results. The caller sets, before including:

     spec_tag     short identifier, e.g. "04 joint / mc_yemp"
     spec_cmd     the command as actually issued, macros already expanded
     spec_rhs     one line describing the regressors
     spec_sample  one line describing the estimation sample
     spec_wt      one line describing the weights
     spec_ref     one line describing the reference / recentring
     spec_out     the file this model's coefficients are written to

   Everything numeric is read from e(), so the log reports what Stata actually
   estimated rather than what the comments claim. The pipeline entry points
   truncate ${MC_OUT}/../audit/REGRESSIONS.txt before stage 04 runs.        */
capture mkdir "${H_ROOT}/outputs/audit"
tempname _ls
file open `_ls' using "${H_ROOT}/outputs/audit/REGRESSIONS.txt", write append text
file write `_ls' "MODEL  `spec_tag'" _n
file write `_ls' "  command    `spec_cmd'" _n
file write `_ls' "  regressors `spec_rhs'" _n
file write `_ls' "  absorbed   county_id (county fixed effects)" _n
file write `_ls' "  sample     `spec_sample'" _n
file write `_ls' "  weights    `spec_wt'" _n
file write `_ls' "  vce        cluster(county_id)" _n
file write `_ls' "  fitted     N=" %12.0fc (e(N)) "  clusters=" %6.0f (e(N_clust)) "  df_r=" %6.0f (e(df_r)) _n
file write `_ls' "  reference  `spec_ref'" _n
file write `_ls' "  writes     `spec_out'" _n
file write `_ls' "" _n
file close `_ls'
