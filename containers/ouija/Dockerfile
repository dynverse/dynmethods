FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_github("kieranrcampbell/ouija")'
RUN R -e 'devtools::install_cran("rstan")'
RUN R -e 'devtools::install_cran("coda")'

ADD . /code

ENTRYPOINT /code/run.sh
