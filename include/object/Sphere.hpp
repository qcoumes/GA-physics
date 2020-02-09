#ifndef GA_PHYSICS_SPHERE_HPP
#define GA_PHYSICS_SPHERE_HPP

#include <vector>

#include <GL/glew.h>
#include <c3ga/Mvec.hpp>

#include <object/ObjectVertex.hpp>


namespace object {
    
    class Sphere {
        
        private:
            static constexpr GLuint VERTEX_ATTR_POSITION = 0;
            static constexpr GLuint VERTEX_ATTR_NORMAL = 1;
            static constexpr GLuint VERTEX_ATTR_TEXTURE = 2;
            
        private:
            std::vector<ObjectVertex> vertices;
            c3ga::Mvec<GLfloat> multivector;
            GLboolean modified;
            GLuint64 longitude;
            GLuint64 latitude;
            GLuint vbo = 0;
            GLuint vao = 0;
            
            void buildVertices();
            
            void fillBuffer();
            
        public:
            
            Sphere(const c3ga::Mvec<GLfloat> &multivector, GLuint64 longitude, GLuint64 latitude);
            
            ~Sphere();
            
            [[nodiscard]] const ObjectVertex *data() const;
            
            [[nodiscard]] GLsizei getVertexCount() const;
            
            void update();
            
            void render() const;
    };
}

#endif // GA_PHYSICS_SPHERE_HPP
