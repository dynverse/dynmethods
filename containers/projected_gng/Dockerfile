FROM dynverse/dynwrap:r

RUN apt-get install -y libudunits2-dev

RUN R -e 'devtools::install_github("rcannood/gng")'

ADD . /code

ENTRYPOINT /code/run.sh
