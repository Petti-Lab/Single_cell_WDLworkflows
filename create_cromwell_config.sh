#!/bin/bash

set -o errexit

function usage
{
    echo ""
    echo "usage: create_cromwell_config -o <output.config> -d <outdir> -l <logdir> -G <compute_group> -g <lsf_group> -v <volumes> -h"
    echo ""
    echo "  -o | --output_file    output file (required)"
    echo "  -d | --outdir         output directory for workflow (required)"
    echo "  -l | --logdir         logfile directory for workflow (required)" 
    echo "  -G | --compute-group  compute group to use for LSF jobs. If not provided,"
    echo "                        will look in the LSF_COMPUTE_GROUP env variable"
    echo "  -q | --compute-queue  LSF queue to run jobs in. If not provided, will"
    echo "                        look in the LSF_COMPUTE_QUEUE env variable"
    echo "  -g | --lsf-job-group  LSF job group to place jobs into" 
    echo "  -v | --volumes        disk volumes to mount in LSF jobs. If not provided,"
    echo "                        will look in the LSF_DOCKER_VOLUMES env variable"
    echo "  -h | --help           show help/options"
    echo ""
    exit 1
}

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    usage
fi

while [[ ! -z "$1" ]]; do
    case $1 in
        -o | --outfile )        shift
                                outfile=$1
                                ;;
        -d | --outdir )         shift
                                outdir=$1
                                ;;
        -l | --logdir )         shift
                                logdir=$1
                                ;;
        -G | --compute-group )  shift
                                compute_group=$1
                                ;;
        -q | --compute-queue )  shift
                                compute_queue=$1
                                ;;
        -g | --lsf-job-group )  shift
                                lsf_group=$1
                                ;;
        -v | --volumes )        shift
                                volumes=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

#resolve full paths, make dirs if needed
if [[ ! -d $outdir ]];then
    mkdir -p $outdir/cromwell-db
fi
outdir=$(readlink -f $outdir)

if [[ ! -d $logdir ]];then
    mkdir -p $logdir
fi
outdir=$(readlink -f $outdir)

#if volumes not provided, then look for env variable
if [[ -z $volumes ]];then
    if [[ ! -z $LSF_DOCKER_VOLUMES ]];then
        volumes="$LSF_DOCKER_VOLUMES"
    fi
fi 
#if compute group not provided, then look for env variable
if [[ $compute_group == "" ]];then
    if [[ ! -z $LSF_COMPUTE_GROUP ]];then
        compute_group=$LSF_COMPUTE_GROUP
    else
        echo "compute group must be provided if LSF_COMPUTE_GROUP env variable is not set"
        exit 1;
    fi
fi
#if compute queue not provided, then look for env variable
if [[ $compute_queue == "" ]];then
    if [[ ! -z $LSF_COMPUTE_QUEUE ]];then
        compute_queue=$LSF_COMPUTE_QUEUE
    else
        echo "compute queue must be provided if LSF_COMPUTE_QUEUE env variable is not set"
        exit 1;
    fi
fi

outdir=$(readlink -f $outdir)
logdir=$(readlink -f $logdir)

echo 'include required(classpath("application"))' >$outfile
echo '' >>$outfile
echo 'backend {' >>$outfile
echo '  default = "LSF"' >>$outfile
echo '  providers {' >>$outfile
echo '    LSF {' >>$outfile
echo '      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"' >>$outfile
echo '      config {' >>$outfile
echo '        runtime-attributes = """' >>$outfile
echo '        Int cpu = 1' >>$outfile
echo '        Int memory_kb = 4096000' >>$outfile
echo '        Int memory_mb = 4096' >>$outfile
echo '        String? docker' >>$outfile
echo '        """' >>$outfile
echo '' >>$outfile
echo '        submit-docker = """' >>$outfile
echo '        LSF_DOCKER_VOLUMES='\''${cwd}:${docker_cwd} '$volumes''\'' \' >>$outfile
echo '        LSF_DOCKER_PRESERVE_ENVIRONMENT=false \' >>$outfile
echo '        bsub \' >>$outfile
echo '        -J ${job_name} \' >>$outfile
echo '        -cwd ${cwd} \' >>$outfile
echo '        -o '$logdir'/cromwell-%J.out \' >>$outfile
echo '        -e '$logdir'/cromwell-%J.err \' >>$outfile
echo '        -q '$compute_queue' \' >>$outfile
echo '        -G '$compute_group' \' >>$outfile
if [[ ! -z "$lsf_job_group" ]];then
    echo '        -g '$lsf_job_group' \' >>$outfile
fi
echo '        -a "docker0(${docker})" \' >>$outfile
echo '        -M ${memory_mb}M \' >>$outfile
echo '        -n ${cpu} \' >>$outfile
echo '        -R "span[hosts=1] select[mem>${memory_mb}M] rusage[mem=${memory_mb}M]" \' >>$outfile
echo '        /bin/bash ${script}' >>$outfile
echo '        """' >>$outfile
echo '' >>$outfile
echo '        kill = "bkill ${job_id}"' >>$outfile
echo '        kill-docker = "bkill ${job_id}"' >>$outfile
echo '        check-alive = "bjobs -noheader -o stat ${job_id} | /bin/grep '\''PEND\\|RUN'\''"' >>$outfile
echo '        job-id-regex = "Job <(\\d+)>.*"' >>$outfile
echo '        root = "'$outdir'"' >>$outfile
echo '        exit-code-timeout-seconds = 600' >>$outfile
echo '        filesysytems {' >>$outfile
echo '          local {' >>$outfile
echo '            caching {' >>$outfile
echo '              duplication-strategy: [' >>$outfile
echo '                "hard-link", "soft-link", "copy"' >>$outfile
echo '              ]' >>$outfile

echo '              hashing-strategy: "xxh64"' >>$outfile
echo '              fingerprint-size: 10485760' >>$outfile
echo '              check-sibling-md5: false' >>$outfile
echo '            }' >>$outfile
echo '          }' >>$outfile
echo '        }' >>$outfile
echo '      }' >>$outfile
echo '    }' >>$outfile
echo '  }' >>$outfile
echo '}' >>$outfile
echo 'workflow-options {' >>$outfile
echo '  workflow-log-dir = "'$logdir'"' >>$outfile
echo '}' >>$outfile
echo 'database {' >>$outfile
echo '  profile = "slick.jdbc.HsqldbProfile$"' >>$outfile
echo '  db {' >>$outfile
echo '    driver = "org.hsqldb.jdbcDriver"' >>$outfile
echo '    url = """' >>$outfile
echo '    jdbc:hsqldb:file:'$outdir'/cromwell-db;' >>$outfile
echo '    shutdown=false;' >>$outfile
echo '    hsqldb.default_table_type=cached;hsqldb.tx=mvcc;' >>$outfile
echo '    hsqldb.result_max_memory_rows=10000;' >>$outfile
echo '    hsqldb.large_data=true;' >>$outfile
echo '    hsqldb.applog=1;' >>$outfile
echo '    hsqldb.lob_compressed=true;' >>$outfile
echo '    hsqldb.script_format=3' >>$outfile
echo '    """' >>$outfile
echo '    connectionTimeout = 120000' >>$outfile
echo '    numThreads = 1' >>$outfile
echo '   }' >>$outfile
echo '}' >>$outfile
echo 'call-caching {' >>$outfile
echo '  enabled = true' >>$outfile
echo '  invalidate-bad-cache-results = true' >>$outfile
echo '}' >>$outfile
