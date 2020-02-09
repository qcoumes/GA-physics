#ifndef GA_PHYSICS_OBJECTVERTEX_HPP
#define GA_PHYSICS_OBJECTVERTEX_HPP

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>


namespace object {
    
    struct ObjectVertex {
        glm::vec3 position;
        glm::vec3 normal;
        glm::vec2 texture;
        
        ObjectVertex() = default;
        
        ObjectVertex(const glm::vec3 &position, const glm::vec3 &normal, const glm::vec2 &texture);
    };
}

#endif //GA_PHYSICS_OBJECTVERTEX_HPP
