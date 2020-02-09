#include <app/object/ObjectVertex.hpp>


namespace app::object {
    
    ObjectVertex::ObjectVertex(const glm::vec3 &t_position, const glm::vec3 &t_normal, const glm::vec2 &t_texture) :
            position(t_position), normal(t_normal), texture(t_texture) {
    }
}
