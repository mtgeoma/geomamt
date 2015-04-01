# -*- mode: Shell-script; -*-

# sysexits
EX_OK=0           # successful termination
EX__BASE=64       # base value for error messages
EX_USAGE=64       # command line usage error
EX_DATAERR=65     # data format error
EX_NOINPUT=66     # cannot open input
EX_NOUSER=67      # addressee unknown
EX_NOHOST=68      # host name unknown
EX_UNAVAILABLE=69 # service unavailable
EX_SOFTWARE=70    # internal software error
EX_OSERR=71       # system error (e.g., can't fork)
EX_OSFILE=72      # critical OS file missing
EX_CANTCREAT=73   # can't create (user) output file
EX_IOERR=74       # input/output error
EX_TEMPFAIL=75    # temp failure; user is invited to retry
EX_PROTOCOL=76    # remote error in protocol
EX_NOPERM=77      # permission denied
EX_CONFIG=78      # configuration error
EX__MAX=78        # maximum listed value

check_file ()
# check if file exist and have reading permission
# $1: file to be checked
{
    file="$1"
    if [[ -r "$file" ]] ; then
	return $EX_OK
    elif [[ -f "$file" ]] ; then
	echo "file \"$file\" without reading permission" >&2
	exit $EX_NOPERM
    else
	echo "file \"$file\" don't exist" >&2
	exit $EX_NOINPUT
    fi
}

needed_commands ()
# check if commands exist in PATH
# $1: commands list to be checked
{
    commands="$1"
    missing_counter=0
    for needed_command in $commands; do
	if ! hash "$needed_command" >/dev/null 2>&1; then
	    printf "Command not found in PATH: %s\n" "$needed_command" >&2
	    ((missing_counter++))
	fi
    done

    if ((missing_counter > 0)); then
	printf "Minimum %d commands are missing in PATH, aborting\n" \
	    "$missing_counter" >&2
	return $EX_OSFILE
    fi
    return $EX_OK
}
