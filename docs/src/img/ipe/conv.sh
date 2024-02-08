#!/bin/bash
for i in *.pdf
do
    pdf2svg ${i%.*}.pdf ../${i%.*}.svg all
done
