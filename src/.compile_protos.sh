#!/bin/bash
PROTO_SRC=$PWD/protos/
for elem in $PROTO_SRC*.proto; do protoc -I$PROTO_SRC --cpp_out=$PWD $elem; done

