#! /bin/sh

for file in `grep -l 'TRoot' *`
do
  echo $file
  sed 's/TRoot/ExRoot/g' ${file} > ${file}.tmp
  mv ${file}.tmp ${file}
done
