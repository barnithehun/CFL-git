#!/bin/bash
# wait_wave.sh <wave> [<n>]: returns when n candidates of the wave have finished (status done or error in DIAG.csv),
# or when no Julia process is left running. Prints the finished count and each finished id with its status.
wave=$1; n=${2:-0}
cd "$(dirname "$0")"
ids=$(awk -F, -v w="$wave" 'NR>1 && $2==w {print $1}' candidates.csv)
total=$(echo "$ids" | wc -w)
[ "$n" -eq 0 ] && n=$total
while true; do
  fin=0
  for id in $ids; do
    f=runs/$id/DIAG.csv
    [ -f "$f" ] && sed -n 2p "$f" | grep -qE "^[^,]*,[^,]*,[^,]*,(done|error)," && fin=$((fin+1))
  done
  [ $fin -ge $n ] && break
  running=$(tasklist //FI "IMAGENAME eq julia.exe" 2>/dev/null | grep -c julia.exe)
  [ "$running" -eq 0 ] && { echo "no julia process left"; break; }
  sleep 60
done
echo "wave $wave: $fin of $total finished"
for id in $ids; do
  f=runs/$id/DIAG.csv
  [ -f "$f" ] && echo "$id $(sed -n 2p "$f" | cut -d, -f4)"
done
