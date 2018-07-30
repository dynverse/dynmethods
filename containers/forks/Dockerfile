FROM library/python

RUN pip install scipy numpy pandas matplotlib
RUN pip install sklearn seaborn Cython
RUN pip install hdbscan

RUN git clone https://github.com/macsharma/FORKS.git

ADD . /code
ENTRYPOINT python /code/run.py
