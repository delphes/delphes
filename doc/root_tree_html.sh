#! /bin/sh

if [ $# -ne 1 ]
then
  echo " Usage: $0 output_file"
  echo " output_file - output file in HTML format."
  exit
fi

awk '
  BEGIN {
    print "<!doctype html>"
    print "<html>"

    print "<head>"
    print "  <meta HTTP-EQUIV=\"Content-Type\" CONTENT=\"text/html; charset=iso-8859-1\">"
    print "  <meta NAME=\"keywords\" CONTENT=\"root, tree, ntuple, format, description\">"
    print "  <title>root tree description</title>"
    print "  <style>"
    print "    body { font-family: sans-serif; max-width: 800px; line-height: 1.4; margin-left: auto; margin-right: auto; padding: 0em 1em 3em 1em; }"
    print "    h1 { font-weight: normal; }"
    print "    td, th { border: 1px solid #d3d3d3; padding: 0.1em 0.3em; }"
    print "    th { text-align: left; background-color: #f5f5f5; }"
    print "    tr.even td { background-color: #f5f5f5; }"
    print "    tr.class td { font-weight: bold; padding-top: 1em; }"
    print "  </style>"
    print "</head>"

    print "<body>"

    print "<h1>ROOT Tree Description</h1>"
    print "<p>Description of all classes used to store output data.</p>"

    print "<table>"
    print "<tr><th>Parameter</th>"
    print "<th>Definition</th>"
    print "<th>How it was calculated</th></tr>"
  }

  function print_line(name, comment, even, end) {
    if(name != ""){
      if(even) print "<tr class=\"even\">"
      else print "<tr class=\"odd\">"
      print "  <td>"name"</td>"
      split(comment, a, "|");
      print "  <td>"a[1]"</td>"
      print "  <td>"a[2]"</td>"
      print "</tr>"
    }
  }

  /^ *class /{
    print_line(name, comment, even, 1);
    even = 0;
    name = "";
    comment = "";
    split($2, a, ":");
    if(a[1] == "Candidate" || a[1] == "DelphesFactory;") next;
    print "<tr class=\"class\"><td colspan=\"3\" id=\""a[1]"\">class "a[1]"</td></tr>"
  }

  /: public /{
    if($4 == "TObject" || $4 == "SortableObject") next;
    name = sprintf("<a href=\"#%s\">%s</a>", $4, $4);
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
    print "</table>"
    print "</body></html>"
  }' `dirname $0`/../classes/DelphesClasses.h > $1
