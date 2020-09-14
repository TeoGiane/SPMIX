#ifndef SRC_RECORDIO_HPP
#define SRC_RECORDIO_HPP

#include <vector>
#include <fstream>
#include <string>
#include <deque>
#include <streambuf>
#include <iostream>
#include <sstream>
#include <google/protobuf/text_format.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/util/delimited_message_util.h>

#include <errno.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>

template <typename T>
bool writeManyToFile(std::deque<T> messages, std::string filename) {
    int outfd = open(filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0777);
    google::protobuf::io::FileOutputStream fout(outfd);
    int cnt = 0;
    bool success;
    for (auto msg: messages) {
        // msg.PrintDebugString();
        success = google::protobuf::util::SerializeDelimitedToZeroCopyStream(
            msg, &fout);
        if (! success) {
            std::cout << "Writing Failed" << std::endl;
            break;
        }
        cnt += 1;
    }
    fout.Close();
    close(outfd);
    std::cout << "Serialized " << cnt << " messages" << std::endl;
    return success;
}

template <typename T>
std::deque<T> readManyFromFile(std::string filename) {
    int infd = open(filename.c_str(), O_RDWR);
    if (infd == -1)
        std::cout << "errno: " << strerror(errno) << std::endl;
    std::cout << "infd: " << infd << std::endl;

    google::protobuf::io::FileInputStream fin(infd);
    bool keep = true;
    bool clean_eof = true;
    std::deque<T> out;

    while (keep) {
        T msg;
        keep = google::protobuf::util::ParseDelimitedFromZeroCopyStream(
            &msg, &fin, nullptr);
        if (keep)
            out.push_back(msg);
    }
    fin.Close();
    close(infd);
    std::cout << "Deserialized " << out.size() << " messages" << std::endl;
    return out;
}

template <typename T>
T loadTextProto(std::string filename) {
    std::ifstream ifs(filename);
    google::protobuf::io::IstreamInputStream iis(&ifs);
    T out;
    auto success = google::protobuf::TextFormat::Parse(&iis, &out);
    if (! success)
        std::cout << "An error occurred in 'loadTextProto'; success: " <<
        success << std::endl; 
    return out;
}


#endif // SRC_RECORDIO_HPP
