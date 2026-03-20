. ~skel/std.profile

# This is needed (and recommended) so that remote ssh sessions have access
# to aliases, functions, and variables defined in ~/.bashrc.
if [ -f ~/.bashrc ]; then
    . ~/.bashrc
fi
