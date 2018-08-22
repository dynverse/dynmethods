FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("monocle")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
