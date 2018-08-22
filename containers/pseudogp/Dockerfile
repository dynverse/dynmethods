FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_github("kieranrcampbell/pseudogp")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
