FROM library/python

# install python dependencies
RUN pip install scipy pandas numpy matplotlib

# install gpy/gpcluster/gpfates
RUN pip install GPy
RUN pip install git+https://github.com/SheffieldML/GPclust.git
RUN pip install git+https://github.com/Teichlab/GPfates.git@bccd5496b4121b3e634ce7cd5b0bff823b2850fa

ADD . /code
ENTRYPOINT python /code/run.py
