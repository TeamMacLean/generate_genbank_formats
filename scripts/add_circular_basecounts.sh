#!/bin/bash
basecount=$2
cat $1 | sed  "s/DNA        /DNA     circular/;s/ORIGIN/BASE COUNT     ${basecount}\nORIGIN/"
