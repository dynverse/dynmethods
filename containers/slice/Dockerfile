FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_github("dynverse/SLICE")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
