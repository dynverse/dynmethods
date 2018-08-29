FROM dynverse/dynwrap:bioc

LABEL version 0.1.0.1

RUN R -e 'devtools::install_cran("sincell")'

ADD . /code

ENTRYPOINT /code/run.sh
