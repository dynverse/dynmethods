FROM library/python
RUN pip install pandas numpy matplotlib
RUN pip install pymatcher

ADD . /code

ENTRYPOINT python /code/run.py
