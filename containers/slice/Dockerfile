FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_github("dynverse/SLICE")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
