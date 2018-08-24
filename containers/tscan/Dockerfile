FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("TSCAN")'

ADD . /code

ENTRYPOINT /code/run.sh
