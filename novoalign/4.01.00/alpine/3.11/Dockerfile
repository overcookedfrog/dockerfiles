FROM alpine:3.11 as builder
ADD novocraftV4.01.00.Linux3.10.0.tar.gz .
RUN apk add --update build-base curl zlib-dev bzip2-dev ncurses-dev xz-dev curl-dev
WORKDIR /app
RUN curl -L https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 | tar xjf - --strip-components=1 \
    && ./configure \
    && make -j \
    && make install

FROM alpine:3.11
RUN apk add --no-cache procps zlib libcurl libbz2 xz-libs ncurses bash
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=builder /novocraft/novoalign   /usr/local/bin/novoalign
COPY --from=builder /novocraft/novosort    /usr/local/bin/novosort  
COPY --from=builder /novocraft/novoindex   /usr/local/bin/novoindex 
COPY --from=builder /novocraft/novoutil    /usr/local/bin/novoutil  
