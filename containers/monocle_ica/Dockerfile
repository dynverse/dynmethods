FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("monocle")'

ADD . /code

ENTRYPOINT /code/run.sh
