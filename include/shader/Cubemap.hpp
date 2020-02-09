#ifndef GA_PHYSICS_CUBEMAP_HPP
#define GA_PHYSICS_CUBEMAP_HPP

#include <memory>

#include <GL/glew.h>

#include <util/INonCopyable.hpp>
#include <util/Image.hpp>


namespace mastercraft::shader {
    
    class Cubemap : public util::INonCopyable {
        private:
            GLuint textureId = 0;
        
        public:
            
            Cubemap() = default;
            
            explicit Cubemap(std::unique_ptr<util::Image> texture[6]);
            
            [[nodiscard]] GLuint getCubemapId() const;
    };
}

#endif // GA_PHYSICS_CUBEMAP_HPP
