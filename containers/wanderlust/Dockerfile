FROM library/python

# install python dependencies
RUN pip install pandas numpy matplotlib
RUN pip install Cython

# install wishbone
RUN pip install git+https://github.com/dynverse/pywishbone --upgrade --upgrade-strategy only-if-needed

ADD . /code
ENTRYPOINT python /code/run.py
