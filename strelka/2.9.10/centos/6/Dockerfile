FROM centos:6
ENV PATH=/app/bin:$PATH
WORKDIR /app
RUN curl -L https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 | tar xjf - --strip-components=1
