#! /bin/sh

if [ $# -ne 1 ]
then
  echo " Usage: $0 config_file"
  echo " config_file - configuration file in Tcl format."
  exit
fi

awk '
  BEGIN{
    print "digraph G {"
    print "rankdir=LR";
  }
  $1~/module/{
    dst=$3
  }
  $2~/InputArray/{
    split($3, src, "/");
    if(src[1] == "Delphes") src[1] = "Reader";
    print dst"[shape=box, style=rounded];";
    print src[1]"->"dst" [label="src[2]"];"
  }
  $2~/Branch/{
    split($3, src, "/");
    if(src[1] == "Delphes") src[1] = "Reader";
    print "\"Branch "$4"\" [shape=box, style=\"rounded,filled\", fillcolor=lightgrey];";
    print src[1]"->\"Branch "$4"\" [label="src[2]"];";    
  }
  END{
    print "}"
  }' $1 | dot -Tpng > data_flow.png
