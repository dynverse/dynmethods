FROM dynverse/dynwrap:bioc

RUN apt-get install -y libgsl-dev

RUN R -e 'install.packages("cellTree")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
