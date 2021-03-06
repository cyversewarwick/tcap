# run as:
#   docker build -t tcap .

# base everything on a recent Ubuntu
FROM debian:latest

# get system packages up to date then install a basic scientific python
RUN apt-get update && apt-get -y upgrade && \
    apt-get -y install python3 \
         python3-numpy python3-scipy python3-pandas python3-matplotlib python3-seaborn \
         ttf-bitstream-vera

# add and configure our code
RUN mkdir /scripts
COPY scripts/ /scripts/

# grab the expression filter utility
RUN apt-get -y install git
RUN git clone https://github.com/cyversewarwick/expression_filter
RUN cp /expression_filter/scripts/expression_filter.py /scripts

MAINTAINER Krzysztof Polanski <k.t.polanski@warwick.ac.uk>

#set up analysis crash text file
RUN apt-get -y install git
RUN git clone https://github.com/cyversewarwick/analysis_crash.git

# so this is what is going to run by default when you trigger this, in the virtual machine
# call the wigwams wrapper from the other directory while staying in /agave with the files
ENTRYPOINT ["bash", "/scripts/tcap_tarwrapper.sh"]