FROM library/python

# install common dependencies
RUN pip install scipy numpy pandas matplotlib

# install common dependencies
RUN pip install sklearn seaborn Cython

RUN pip install h5py flask gunicorn numexpr scikit-learn simplejson psutil

RUN pip install pot

RUN pip install wot

RUN pip install --upgrade pip

# upgrade to latest version of wot
RUN pip install --upgrade git+https://github.com/broadinstitute/wot.git@develop#egg=wot

ADD . /code

ENTRYPOINT python /code/run.py
