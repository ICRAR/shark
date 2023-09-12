reference_model_name=$1
shift
this_model_name=$1
shift

python scripts/compare_galaxies.py -m \
  "${TEST_OUTPUTS_DIR}/${TEST_SIM_NAME}/${reference_model_name}/199/0/galaxies.hdf5" \
  "${TEST_OUTPUTS_DIR}/${TEST_SIM_NAME}/${this_model_name}/199/0/galaxies.hdf5" \
  $@
