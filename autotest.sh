#!/bin/bash

function isOk {
  if [ $1 == 0 ] ; then echo OK ; else echo ERROR ; fi
  }

echo -n "test 1... "
./bipype -m run -r /home/pszczesny/hajduk/test_samples --out_dir /home/pszczesny/hajduk/bipype/output > /dev/null 2> /dev/null
isOk $?
echo -n "test 2... " 
./bipype >/dev/null 2>/dev/null
isOk $?
