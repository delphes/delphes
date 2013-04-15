#! /bin/sh

awk '
  BEGIN {  
    print "||= Parameter =||= Definition =||= How it was calculated =||";
  }

  function print_line(name, comment, even, end) {
    if(name != ""){
      split(comment, a, "|");
      print "|| "name" || "a[1]" || "a[2]" ||";
    }
  }

  /^ *class /{
    print_line(name, comment, even, 1);
    even = 0;
    name = "";
    comment = "";
    split($2, a, ":");
    if(a[1] == "Candidate" || a[1] == "DelphesFactory;") next;
    print "|| || || ||";
    print "|||| class "a[1]" || ||"
    print "|| || || ||"
  }

  /: public [^S]/{
    if($4 == "TObject") next;
    name = $4;
    split($2, a, ":");
    comment = sprintf("%s inherits all %s parameters", a[1], $4);
  }

  /^ *[A-Za-z_]* [A-Za-z].*; \/\/ / {
    print_line(name, comment, even, 0);
    split($2, a, ";");
    name = a[1];
    split($0, a, "// ");
    comment = a[2];
    even = !even;
  }

  /^ +\/\/ /{split($0, a, "// "); comment = comment" "a[2]}
  END {
    print_line(name, comment, even, 1);
  }' ../classes/DelphesClasses.h

