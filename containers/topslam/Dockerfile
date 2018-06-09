FROM library/python:2.7

# install python dependencies
RUN pip install pandas numpy matplotlib

# install GPy
RUN git clone https://github.com/SheffieldML/GPy.git ; \
    cd GPy; \
    python setup.py develop

# install GPy
RUN git clone https://github.com/mzwiessele/topslam.git ; \
    cd topslam; \
    python setup.py develop

ADD . /code
ENTRYPOINT python /code/run.py
