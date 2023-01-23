#!/usr/bin/env bash
######################################################################
# test runner
######################################################################

# (0) paths
base=$(git rev-parse --show-toplevel)
if [ ! -d ${base} ]; then
    echo Cannot find the base directory.
    exit 1
fi
data=${base}/dat
inst=${base}/instances
test=${base}/test

source ${base}/src/scripts/functions.sh
alloptions="$*"

# options for parallel
popt=""
dry_run=false
instances="wapI"

function usage() {
    echo "$(basename $0) [-n] [-j n] [-t timelimit] [-i instances] options*"
    exit 1
}

## (1) options
## https://gist.github.com/magnetikonline/0e44ab972a7efa3ac138
while getopts ":nhj:i:t:" o; do
    case ${o} in
        n) popt="${popt} --dry-run"; dry_run=true;;
        j) popt="${popt} -j ${OPTARG}";;
        t) timelimit="${OPTARG}";;
	i) instances=${OPTARG};;
	h) usage; exit 0;;
    esac
done
shift $((OPTIND-1))

if [ $# -lt 2 ]; then
    echo Usage: $(basename $0) [options] test-number description [additional options]
    exit 1
fi
testno=$(printf %03d $1)
outdir=${test}/test${testno}
desc=$2
shift 2
options=$*

# (3) prepare output directory
prepare_out ${outdir} ipsolve
solver=${outdir}/ipsolve
cd ${outdir}

h_options=""
if [ ! -z ${timelimit+x} ]; then
    h_options="--timelimit ${timelimit}"
fi

function run_setI() {
    parallel ${popt} "${solver} ${h_options} --instance {1} ${options}" ::: ${data}/Battarra,etal/*.txt
}

function run_setII() {
    parallel ${popt} "${solver} ${h_options} --instance {1} ${options}" ::: ${data}/setII/*.txt
}

function run_setIII() {
    parallel ${popt} "${solver} ${h_options} --instance {1} ${options}" ::: ${data}/setIII/*.txt
}

function run_ALWABP() {
    parallel ${popt} "${solver} ${h_options} --instance ${data}/alwabp/alwabp/{1}/{2} ${options}" ::: heskia roszieg tonge wee-mag ::: {01..80}
}

function run_set() {
    local instances="$1"

    (
	if [ "${instances}" = "wapI" ]; then
	    run_setI
	elif  [ "${instances}" = "wapII" ]; then
	    run_setII
	elif  [ "${instances}" = "wapIII" ]; then
	    run_setIII
	else ## ALWABP
	    run_ALWABP
	fi
    ) > ${outdir}/solver.raw
    numtags=$(awk '{print $1}' solver.raw | sort | uniq | wc -l)
    if [ $numtags -eq 1 ]; then
	column -t ${outdir}/solver.raw > ${outdir}/solver.log
    elif [ $numtags -eq 2 ]; then
	paste - - < ${outdir}/solver.raw | column -t > ${outdir}/solver.log
    elif [ $numtags -eq 3 ]; then
	paste - - - < ${outdir}/solver.raw | column -t > ${outdir}/solver.log
    elif [ $numtags -eq 4 ]; then
	paste - - - - < ${outdir}/solver.raw | column -t > ${outdir}/solver.log
    else
	paste - - - - - < ${outdir}/solver.raw | column -t > ${outdir}/solver.log
    fi
}

# (4) run
save_config ${outdir}
run_set ${instances}
finish ${outdir}
