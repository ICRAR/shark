model_name=$1
shift

${CI_BUILD_DIR}/shark ${TEST_CONFIG_FILE} \
    -o execution.output_bh_histories=true -o execution.snapshots_bh_histories=199 \
    -o simulation.redshift_file=${TEST_INPUTS_DIR}/redshifts.txt \
    -o simulation.tree_files_prefix=${TEST_INPUTS_DIR}/tree_199 \
    -o simulation.sim_name=${TEST_SIM_NAME} \
    -o execution.name_model=$model_name \
    $@
