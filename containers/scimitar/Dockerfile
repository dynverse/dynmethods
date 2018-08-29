FROM dynverse/dynwrap:py2.7

LABEL version 0.1.0.1

# install R before rpy2
RUN apt-get update && apt-get install -y r-base

RUN pip install rpy2==2.8
RUN pip install munkres
RUN pip install git+https://github.com/dimenwarper/pyroconductor
RUN pip install git+https://github.com/dimenwarper/scimitar

# install R dependencies
RUN R -e "install.packages('corpcor', repos = 'http://cran.us.r-project.org')"

ADD . /code
ENTRYPOINT /code/run.sh
