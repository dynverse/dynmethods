FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_github("kieranrcampbell/mfa")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
