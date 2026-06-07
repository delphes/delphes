import ast

import tkinter as tk

from argparse import ArgumentParser
from itertools import batched, groupby
from pathlib import Path

script = r"""
proc format_string {s} {
    set lines {}
    foreach l [split $s \n] {
        set l [string trim $l]
        if {$l ne ""} {
            lappend lines $l
        }
    }
    set n [llength $lines]
    if {$n > 1} {
        set tmp [join $lines \n]
        return '''$tmp'''
    } elseif {$n == 1} {
        set tmp [lindex $lines 0]
        set tmp [regsub -all {\s} $tmp {}]
        return '$tmp'
    } else {
        return ''
    }
}

proc format_dict {value} {
    set temp_vars [list EtaBins PhiBins eta i pi]
    regexp {^value is a (.*?) with a refcount} [::tcl::unsupported::representation $value] -> type
    switch $type {
        dict {
            set tmp {}
            foreach {k v} [dict map {k v} $value {format_dict $v}] {
                if {$k ni $temp_vars} {
                    lappend tmp '$k':$v
                }
            }
            return \{[join $tmp ,]\}
        }
        list {
            set tmp [lmap v $value {format_dict $v}]
            return \[[join $tmp ,]\]
        }
        int - double {
            return [expr {$value}]
        }
        booleanString {
            return [expr {$value ? True : False}]
        }
        default {
            if {[string is integer -strict $value]} {
                return [expr {$value}]
            } elseif {[string is double -strict $value]} {
                return [expr {$value}]
            } elseif {[string is boolean -strict $value]} {
                return [expr {$value ? True : False}]
            }
            return [format_string $value]
        }
    }
}

proc module {class name body} {
    namespace eval ::${name} $body
    set ::${name}::Class $class
}

interp alias {} add {} lappend

proc tcl2dict {} {
    global ExecutionPath

    set card_dict [dict create ExecutionPath $ExecutionPath]

    foreach module $ExecutionPath {
        set vars [info vars ${module}::*]

        set names [lmap v $vars {regsub ^::${module}:: $v {}}]

        set module_dict {}

        foreach v $vars n $names {
            dict set module_dict $n [set $v]
        }

        dict set card_dict $module $module_dict
    }
    return [format_dict $card_dict]
}
"""


def optimize_eta_phi_bins(card):
    for module in card.values():
        if not isinstance(module, dict) or "EtaPhiBins" not in module:
            continue
        new_bins = []
        for phi, group in groupby(batched(module["EtaPhiBins"], 2), key=lambda p: p[1]):
            eta_bins = [item for eta, _ in group for item in (eta if isinstance(eta, list) else [eta])]
            new_bins.extend([eta_bins, phi])
        module["EtaPhiBins"] = new_bins


def tcl2dict(card_file):
    def source(source_file):
        source_path = card_dir / source_file
        source_content = source_path.read_text(encoding="utf-8").replace("\\", "")
        tcl.call("safe_interp", "eval", source_content)

    card_dir = card_file.parent
    card_content = card_file.read_text(encoding="utf-8").replace("\\", "")

    tcl = tk.Tcl()

    tcl.createcommand("safe_source", source)

    tcl.eval("interp create -safe safe_interp")
    tcl.eval("interp alias safe_interp source {} safe_source")

    tcl.call("safe_interp", "eval", script)
    tcl.call("safe_interp", "eval", card_content)

    result = tcl.eval("safe_interp eval tcl2dict")

    card = ast.literal_eval(result)

    optimize_eta_phi_bins(card)

    return card


def main():
    parser = ArgumentParser()
    parser.add_argument("card", help="detector card file name")

    args = parser.parse_args()

    card_file = Path(args.card)
    card_dict = tcl2dict(card_file)

    print(f"card_dict = {card_dict}")


if __name__ == "__main__":
    main()
