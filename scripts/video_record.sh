#!/bin/sh
recordmydesktop -windowid `xwininfo |grep "Window id:"|sed -e "s/xwininfo\:\ Window id:\ // ;s/\ .*//"`
