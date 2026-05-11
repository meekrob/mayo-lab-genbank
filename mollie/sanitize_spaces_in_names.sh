#!/usr/bin/env bash

find . -depth -name '* *' -execdir bash -c 'mv -v "$0" "${0// /_}"' {} \;
