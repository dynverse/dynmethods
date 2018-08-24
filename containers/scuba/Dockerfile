FROM dynverse/dynwrap:py3.6

# install R before rpy2
RUN apt-get update
RUN apt-get install -y r-base

# install scuba
RUN pip install rpy2
RUN pip install jinja2
RUN pip install git+https://github.com/dynverse/PySCUBA.git # first install without upgrade of dependencies
RUN pip install git+https://github.com/dynverse/PySCUBA.git --upgrade --no-dependencies --no-cache-dir -U # now upgrade PySCUBA

# install legacy princurve for now
RUN R -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'
RUN R -e 'devtools::install_github("dynverse/princurve@69b85ad4709b15e5b40f8541f4b5e2ca9059be3a")'

ADD . /code

ENTRYPOINT /code/run.sh
