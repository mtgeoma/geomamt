#!/bin/bash
# This script scans the output of "ps -e" every $1 seconds, until grep
# does not find $2 in this output anymore, then bell.sh is called.
#
# Copyright (C) 2004,2005 Yves Gablin (gablin@fr.fm)
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License version 2 as published
# by the Free Software Foundation: http://www.gnu.org/licenses/gpl.html
#
# Parameters:
# $1 : "sleep" parameter (see "sleep" man page)
# $2 : "grep" parameters (this is a regular expression)
#
# Usage examples:
# watch.sh 1m make  # scan every minute
# watch.sh 15 make  # scan every 15 seconds

this="$(basename "$0")"
sleeptime=$1
shift

while [ -n "$(ps -e -o cmd | grep "$@" | grep -vF "$this" | grep -vF grep)" ]; do
        sleep $sleeptime
done
