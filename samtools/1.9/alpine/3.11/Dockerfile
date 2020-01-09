FROM alpine:3.11 as builder
RUN apk add --no-cache build-base curl curl-dev zlib-dev bzip2-dev xz-dev ncurses-dev
WORKDIR /app
RUN curl -L https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 | tar xj --strip-components=1 \
    && ./configure \
    && make -j \
    && make install

FROM alpine:3.11
RUN apk add --no-cache procps libstdc++ zlib libcurl libbz2 xz-libs ncurses bash
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
