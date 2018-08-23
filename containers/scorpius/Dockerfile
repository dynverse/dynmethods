FROM dynverse/dynwrap:r

RUN R -e 'install.packages("SCORPIUS")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
