FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("TSCAN")'

ADD . /code

ENTRYPOINT /code/run.sh
