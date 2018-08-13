FROM dynverse/dynwrap:r

RUN R -e 'install.packages("SLICER")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
