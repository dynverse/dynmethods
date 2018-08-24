FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_github("kieranrcampbell/ouija")'
RUN R -e 'install.packages("rstan")'
RUN R -e 'install.packages("coda")'

ADD . /code

ENTRYPOINT /code/run.sh
