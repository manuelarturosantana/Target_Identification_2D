function comparison = run_all_identification_methods()
%RUN_ALL_IDENTIFICATION_METHODS  Run every numbered method and compare results.

run_identification_00_baseline
run_identification_01_hard_visibility_filter
run_identification_02_residue_comparison

comparison = read_identification_results();
end
