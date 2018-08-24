FROM dynverse/dynwrap:r

RUN R -e 'install.packages("SCORPIUS")'

ADD . /code

ENTRYPOINT /code/run.sh
