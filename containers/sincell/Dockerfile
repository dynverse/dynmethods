FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("sincell")'

ADD . /code

ENTRYPOINT /code/run.sh
