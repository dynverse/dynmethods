FROM dynverse/dynwrap:r

RUN R -e 'devtools::install_github("kieranrcampbell/mfa")'

ADD . /code

ENTRYPOINT /code/run.sh
