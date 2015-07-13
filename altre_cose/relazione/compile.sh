#!/bin/bash

pdflatex main.tex
pdflatex main-frn.tex
pdflatex main.tex

rm main.log main.toc main.lof main.lot main.bcf main.aux
