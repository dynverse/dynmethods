FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("destiny")'

ADD . /code

ENTRYPOINT /code/run.sh
