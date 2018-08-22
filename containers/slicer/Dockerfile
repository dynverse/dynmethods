FROM dynverse/dynwrap:r

RUN R -e 'install.packages("SLICER")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
