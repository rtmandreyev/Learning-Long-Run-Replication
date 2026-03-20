#!/bin/sh

case `uname` in
    Darwin)
        exec /opt/X11/bin/quartz-wm ;;
    Linux)
        . /etc/X11/xinit/xinitrc ;;
esac
