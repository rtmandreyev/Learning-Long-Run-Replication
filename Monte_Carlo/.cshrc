#	source the standard .cshrc file
source	~skel/std.cshrc

# Operating system customizations
set _uname = `uname -s`
if ( $_uname =~ SunOS ) then
	# Set various Solaris variables
	;
else if ($_uname =~ Linux) then
	# Set various Linux variables
	;
else if ($_uname =~ Darwin) then
	# Set various Mac variables
	;
endif
unset _uname
