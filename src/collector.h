#ifndef SRC_COLLECTOR_HPP
#define SRC_COLLECTOR_HPP

#include <deque>
#include <list>
#include <fstream>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <cstdio>
#include <cstring>
#include <cerrno>

#include "recordio.h"

using namespace google::protobuf::io;

template<typename T>
class Collector {
 protected:
     bool in_memory;
     std::string filename;
     std::deque<T> chains;
     std::ofstream fout;

 public:
     Collector(std::string filename): filename(filename) {
         in_memory = false;
     }

     Collector(int max_num_iter) {
         in_memory = true;
         // chains.reserve(max_num_iter);
     }

     ~Collector() {}

     void collect(T state) {
         if (in_memory) {
             chains.push_back(state);
         } else {
             std::string s;
             fout.open(filename, std::ios::out | std::ios::app);
             if (fout.fail())
                throw std::ios_base::failure(std::strerror(errno));

            fout.exceptions(
                fout.exceptions() | std::ios::failbit | std::ifstream::badbit);

            state.SerializeToString(&s);
            fout << s << std::endl;
         }
     }

     T get(int iternum) {
         return chains[iternum];
     }

     void saveToFile(std::string fname) {
         writeManyToFile(chains, filename);
     }

     void loadFromFile(std::string filename) {
         chains = readManyFromFile<T>(filename);
     }
};

#endif
