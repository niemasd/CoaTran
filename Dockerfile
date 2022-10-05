# Minimal Docker image for CoaTran using Alpine base
FROM alpine:3.13.5
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# install CoaTran
RUN apk update && \
    apk add bash g++ make && \
    wget -qO- "https://github.com/niemasd/CoaTran/archive/refs/tags/0.0.3.tar.gz" | tar -zx && \
    cd CoaTran-* && \
    make && \
    mv coatran_* /usr/local/bin/ && \
    cd .. && \
    rm -rf CoaTran-*
