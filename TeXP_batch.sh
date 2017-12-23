#!/bin/bash
INSTALL_DIR="/src/texp"

usage() { printf "\nUsage: $0 [-f input.fastq] [-t <interger>] [-o <output_path>] [-n <string>]\n\n -f: Input file (fastq,fastq.gz,sra)\n -t: Number of threads\n -o: Output path (i.e. ./ or ./processed)\n -n: Sample name (i.e. SAMPLE01)\n" 1>&2; exit 1; }

while getopts ":f:o:t:n:" o; do
    case "${o}" in
        f )
            input_file=${OPTARG}
            ;;
        o )
            output_dir=${OPTARG}
            ;;
        t )
            threads=${OPTARG}
            ;;
		n )
            output_name=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${threads}" ]; then
    threads=1
fi

if [ ! -f $INSTALL_DIR/Makefile ]; then
	echo "ERR: Could not find TeXP Makefile at: $INSTALL_DIR/Makefile please fix set the INSTALL_DIR variable in $0" 1>&2; exit 1;
fi

if [ ! -f $input_file ]; then
	echo "ERR: Could not find INPUT file: \"$input_file\" use -f to properly set the path to the input file" 1>&2; usage ;
fi


if [ -z "${input_file}" ] || [ -z "${output_dir}" ] || [ -z "{threads}" ] || [ -z "{output_name}" ]; then
    usage
fi

tar xvf $input_file -C ./ > packedfiles;
for r_input_file in $(cat packedfiles); do
	echo make -f $INSTALL_DIR/Makefile INPUT_FILE_PATH=$r_input_file OUTPUT_DIR=$output_dir N_THREADS=$threads SAMPLE_NAME=$output_name;
done

