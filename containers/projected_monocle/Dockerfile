FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("monocle")'

ADD . /code

ENTRYPOINT /code/run.sh
