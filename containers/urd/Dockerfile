FROM dynverse/dynwrap:bioc

RUN apt-get install -y libudunits2-dev

RUN R -e 'devtools::install_github("farrellja/URD")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
