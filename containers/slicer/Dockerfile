FROM dynverse/dynwrap:r

RUN R -e 'install.packages("SLICER")'

ADD . /code

ENTRYPOINT /code/run.sh
