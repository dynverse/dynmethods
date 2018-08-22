FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("destiny")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
