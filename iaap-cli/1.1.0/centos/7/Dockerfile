FROM centos:7
ENV PATH=/app/iaap-cli:$PATH
WORKDIR /app
RUN curl -L http://webdata.illumina.com.s3-website-us-east-1.amazonaws.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0-sha.80d7e5b3d9c1fdfc2e99b472a90652fd3848bbc7.tar.gz | tar zxf - --strip-components=1
RUN yum install -y icu
RUN chmod a+x /app/iaap-cli/iaap-cli /app/iaap-cli/createdump
