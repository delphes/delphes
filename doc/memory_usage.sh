#! /bin/sh

exec "$@" &
pid=$!
echo $pid 
while true
do
  ps -p $pid -o vsize= -o rss= >> $pid.txt
  sleep 0.1
done
