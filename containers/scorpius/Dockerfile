FROM dynverse/dynwrap:r

RUN R -e 'install.packages("SCORPIUS")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
