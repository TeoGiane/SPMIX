#!/bin/bash
mkdir -p cpp_proto
PROTO_SRC=../inst/proto
for elem in $PROTO_SRC/*.proto; do protoc -I$PROTO_SRC --cpp_out=$PWD/cpp_proto $elem; done

