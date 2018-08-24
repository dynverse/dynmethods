FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_github("kieranrcampbell/phenopath")'

ADD . /code

ENTRYPOINT /code/run.sh
