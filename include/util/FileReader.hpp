#ifndef GA_PHYSICS_FILEREADER_HPP
#define GA_PHYSICS_FILEREADER_HPP

#include <string>

namespace mastercraft::util {
    
    class FileReader {
        
        public:
            
            FileReader() = delete;
            
            static std::string read(const std::string &path);
    };
}

#endif // GA_PHYSICS_FILEREADER_HPP
