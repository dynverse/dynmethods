FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_cran("SCORPIUS")'

ADD . /code

ENTRYPOINT /code/run.sh
