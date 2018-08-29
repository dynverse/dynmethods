FROM dynverse/dynwrap:r

LABEL version 0.1.0.1

RUN R -e 'devtools::install_github("kieranrcampbell/ouija")'
RUN R -e 'devtools::install_cran("rstan")'
RUN R -e 'devtools::install_cran("coda")'

ADD . /code

ENTRYPOINT /code/run.sh
