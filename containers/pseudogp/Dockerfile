FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_github("kieranrcampbell/pseudogp")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
