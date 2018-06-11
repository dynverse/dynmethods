FROM library/python

# install python dependencies
RUN pip install pandas numpy matplotlib

# install R before rpy2
RUN apt-get update
RUN apt-get install -y r-base

# install scuba
RUN pip install rpy2
RUN pip install jinja2
RUN pip install git+https://github.com/dynverse/PySCUBA.git # first install without upgrade of dependencies
RUN pip install git+https://github.com/dynverse/PySCUBA.git --upgrade --no-dependencies --no-cache-dir -U # now upgrade PySCUBA

# install R dependencies
RUN R -e "install.packages('princurve', repos = 'http://cran.us.r-project.org')"

ADD . /code
ENTRYPOINT python /code/run.py
