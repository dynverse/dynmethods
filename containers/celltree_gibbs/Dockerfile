FROM dynverse/dynwrap:bioc

RUN apt-get install -y libgsl-dev

RUN R -e 'devtools::install_cran("cellTree")'

ADD . /code

ENTRYPOINT /code/run.sh
