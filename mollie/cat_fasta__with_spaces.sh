#!/usr/bin/env bash
find . -name '*.fasta' -exec cat {} + > allseqs.fasta
