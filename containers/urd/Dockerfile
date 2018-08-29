FROM dynverse/dynwrap:bioc

LABEL version 0.1.0.1

RUN apt-get install -y libudunits2-dev

RUN R -e 'devtools::install_github("farrellja/URD")'

ADD . /code

ENTRYPOINT /code/run.sh
