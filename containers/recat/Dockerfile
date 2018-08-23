FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_github("dynverse/reCAT")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
