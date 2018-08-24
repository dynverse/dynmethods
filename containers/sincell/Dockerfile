FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("sincell")'

ADD . /code

ENTRYPOINT /code/run.sh
