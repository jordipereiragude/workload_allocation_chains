######################################################################
# common functions for shell scripts
######################################################################

## globals: desc, popt, alloptions, testno
function save_config() {
    local outdir="$1"
    
    config=${outdir}/config.dat

    echo "Start: $(date)" > ${config}
    echo "Status: Ok" >> ${config}
    echo "Description: ${desc}" >> ${config}
    echo "Parallel options: ${popt}" >> ${config}
    echo "Options: ${alloptions}" >> ${config}
    echo "Test number: ${testno} " >> ${config}
    echo -n "GIT revision: " >> ${config}; git rev-parse HEAD >> ${config}
    echo "hostname: $(hostname)" >> ${config}
    echo "Process limits:" >> ${config}; ulimit -a >> ${config}

    git diff ${base} > ${outdir}/version.diff
}

function prepare_out() {
    local outdir="$1"
    shift

    if [ -d ${outdir} ]; then
	echo ${outdir} exists.
	exit 1
    fi
    mkdir ${outdir}
    
    for b in $*; do
	cp ${base}/src/${b} ${outdir}
    done
}

function finish() {
    local outdir="$1"
    echo "Stop: $(date)" >> ${config}
}
