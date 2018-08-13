FROM dynverse/dynwrap:r

RUN R -e 'install.packages("princurve")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
