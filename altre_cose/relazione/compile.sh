#!/bin/bash

pdflatex main.tex
pdflatex -frn main.tex
pdflatex main.tex

rm main.log main.toc main.lof main.lot main.bcf main.aux
