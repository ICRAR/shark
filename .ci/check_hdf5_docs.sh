# Generate the HDF5 output documentation and check it's up to date
# otherwise tell the user how to update it
hdf5_file="$1"
rst_file="$2"
hdf5_file_basename=`basename "${hdf5_file}"`
full_rst_file="doc/hdf5_properties/${rst_file}"

scripts/properties_as_list.sh "${TEST_OUTPUTS_DIR}/${TEST_SIM_NAME}/${hdf5_file}" > props.rst
_diff="`diff -Naur "${full_rst_file}" props.rst`"
if [ -n "${_diff}" ]; then
	echo "\nThe file ${full_rst_file} is out of date. This probably means that you added a new\n" \
	"dataset to shark's output, but forgot to update the corresponding documentation.\n" \
	"The full difference follows:\n\n${_diff}\n\n" \
	"Please run the script/properties_as_lish.sh script against a ${hdf5_file_basename} file\n" \
	"to re-generate its documentation, then commit your changes. For example:\n\n" \
	"scripts/properties_as_list.sh my-output/model/199/0/${hdf5_file_basename} > ${full_rst_file}" >&2
	exit 1
fi
