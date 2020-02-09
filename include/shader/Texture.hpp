#ifndef GA_PHYSICS_TEXTURE_HPP
#define GA_PHYSICS_TEXTURE_HPP

#include <memory>

#include <GL/glew.h>

#include <util/INonCopyable.hpp>
#include <util/Image.hpp>


namespace mastercraft::shader {
    
    class Texture : public util::INonCopyable {
        private:
            GLuint textureId = 0;
        
        public:
            
            Texture() = default;
            
            explicit Texture(const util::Image *texture);
            
            [[nodiscard]] GLuint getTextureId() const;
    };
}

#endif //GA_PHYSICS_TEXTURE_HPP
