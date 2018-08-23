FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("sincell")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
