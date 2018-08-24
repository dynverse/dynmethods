FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("destiny")'

ADD . /code

ENTRYPOINT /code/run.sh
