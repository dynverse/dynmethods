FROM dynverse/dynwrap:bioc

RUN apt-get install -y libgsl-dev

RUN R -e 'install.packages("cellTree")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
