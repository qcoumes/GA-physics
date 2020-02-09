#include <iostream>
#include <cmath>
#include <vector>

#include <app/object/Sphere.hpp>


namespace app::object {
    
    Sphere::Sphere(const c3ga::Mvec<GLfloat> &t_multivector, GLuint64 t_longitude, GLuint64 t_latitude) :
        multivector(t_multivector), modified(true), longitude(t_longitude), latitude(t_latitude) {
        glGenBuffers(1, &this->vbo);
        glGenVertexArrays(1, &this->vao);
    }
    
    
    Sphere::~Sphere() {
        glDeleteBuffers(1, &this->vbo);
        glDeleteVertexArrays(1, &this->vao);
    }
    
    
    void Sphere::buildVertices() {
        std::vector<ObjectVertex> data;
        glm::vec3 position, normal;
        glm::vec2 texture;
        
        c3ga::Mvec<GLfloat> dual = multivector.dual();
        dual /= dual[c3ga::E0];
        
        GLfloat tmp, c, s, radius = std::sqrt(dual | dual);
        GLfloat rcpLat = 1.f / latitude, rcpLong = 1.f / longitude;
        GLfloat dPhi = 2 * M_PIf32 * rcpLat, dTheta = M_PIf32 * rcpLong;
        
        for (GLuint64 j = 0; j <= longitude; ++j) {
            tmp = j * dTheta;
            c = std::cos(-M_PI_2f32 + tmp);
            s = std::sin(-M_PI_2f32 + tmp);
            
            for (GLuint64 i = 0; i <= latitude; ++i) {
                tmp = i * dPhi;
                
                texture = { i * rcpLat, 1.f - j * rcpLong };
                normal = { std::sin(tmp) * c, s, std::cos(tmp) * c };
                position = radius * normal;
                
                data.emplace_back(position, normal, texture);
            }
        }
        
        for (GLuint64 j = 0; j < longitude; ++j) {
            GLuint64 offset = j * (latitude + 1);
            
            for (GLuint64 i = 0; i < latitude; ++i) {
                this->vertices.push_back(data[offset + i]);
                this->vertices.push_back(data[offset + (i + 1)]);
                this->vertices.push_back(data[offset + latitude + 1 + (i + 1)]);
                this->vertices.push_back(data[offset + i]);
                this->vertices.push_back(data[offset + latitude + 1 + (i + 1)]);
                this->vertices.push_back(data[offset + i + latitude + 1]);
            }
        }
    }
    
    
    void Sphere::fillBuffer() {
        glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(ObjectVertex) * this->vertices.size(), this->vertices.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    
        // Set the VAO
        glBindVertexArray(this->vao);
        glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
        glEnableVertexAttribArray(VERTEX_ATTR_POSITION);
        glEnableVertexAttribArray(VERTEX_ATTR_NORMAL);
        glEnableVertexAttribArray(VERTEX_ATTR_TEXTURE);
        glVertexAttribPointer(
            VERTEX_ATTR_POSITION, 3, GL_FLOAT, GL_FALSE, sizeof(ObjectVertex),
            reinterpret_cast<const GLvoid *>(offsetof(ObjectVertex, position))
        );
        glVertexAttribPointer(
            VERTEX_ATTR_NORMAL, 3, GL_FLOAT, GL_FALSE, sizeof(ObjectVertex),
            reinterpret_cast<const GLvoid *>(offsetof(ObjectVertex, normal))
        );
        glVertexAttribPointer(
            VERTEX_ATTR_TEXTURE, 2, GL_FLOAT, GL_FALSE, sizeof(ObjectVertex),
            reinterpret_cast<const GLvoid *>(offsetof(ObjectVertex, texture))
        );
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }
    
    
    const ObjectVertex *Sphere::data() const {
        return this->vertices.data();
    }
    
    
    GLsizei Sphere::getVertexCount() const {
        return static_cast<GLsizei>(this->vertices.size());
    }
    
    
    void Sphere::update() {
        if (!this->modified) {
            return;
        }
        
        buildVertices();
        fillBuffer();
        
        this->modified = false;
    }
    
    
    void Sphere::render() const {
    }
}
