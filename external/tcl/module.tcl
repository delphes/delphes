
proc module {class name body} {
	eval [format {
		catch {%s} msg opts
		foreach v [info locals] {
			if {$v in {class name body msg opts}} continue
			set ::%s::$v [set $v]
		}
		set ::%s::Class %s
		return {*}$opts $msg
	} $body $name $name $class]
}
