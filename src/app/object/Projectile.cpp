#include <iostream>
#include <cmath>
#include <vector>

#include <c3ga/Tools.hpp>

#include <app/object/Projectile.hpp>
#include <glm/glm.hpp>
#include <app/Config.hpp>
#include <app/Engine.hpp>


#ifndef M_PIf32
#define M_PIf32    3.141592653589793238462643383279502884f /**< PI in floating point precision. */
#define M_PI_2f32  M_PIf32 / 2.f /**< PI/2 in floating point precision. */
#endif

namespace app::object {
    
    Projectile::Projectile(std::shared_ptr<shader::ShaderTexture> shader, std::shared_ptr<shader::Texture> texture) :
            shader(std::move(shader)), texture(std::move(texture)), longitude(20), latitude(20) {
        
        Config *config = Config::getInstance();
        Engine *engine = Engine::getInstance();
        glm::vec3 direction = engine->camera->getFrontVector() * config->getProjSpeed();
        glm::vec3 position = engine->camera->getPosition() + (direction * config->getProjSize());
        
        c3ga::Mvec<GLfloat> a = c3ga::point<GLfloat>(0.5f, 0, 0);
        c3ga::Mvec<GLfloat> b = c3ga::point<GLfloat>(-0.5f, 0, 0);
        c3ga::Mvec<GLfloat> c = c3ga::point<GLfloat>(0, 0.5f, 0);
        c3ga::Mvec<GLfloat> d = c3ga::point<GLfloat>(0, 0, 0.5f);
        
        this->multivector = a ^ b ^ c ^ d;
        this->multivector = fVersor::dilator(config->getProjSize())(this->multivector);
        this->multivector = fVersor::translator(position.x, position.y, position.z)(this->multivector);
        this->bounce = config->getProjBounce();
        this->translator = fVersor::translator(direction.x, direction.y, direction.z);
        
        glGenBuffers(1, &this->vbo);
        glGenVertexArrays(1, &this->vao);
    }
    
    
    Projectile::~Projectile() {
        glDeleteBuffers(1, &this->vbo);
        glDeleteVertexArrays(1, &this->vao);
    }
    
    
    void Projectile::transform(const Versor<GLfloat> &versor) {
        this->multivector = versor(this->multivector);
    }
    
    
    c3ga::Mvec<GLfloat> Projectile::collide(c3ga::Mvec<GLfloat> mvec) {
        return c3ga::Mvec<GLfloat>();
    }
    
    
    bool Projectile::checkCollision() {
        Engine *engine = Engine::getInstance();
        
        c3ga::Mvec<GLfloat> empty = c3ga::Mvec<GLfloat>(), collision = engine->collide(this->multivector);
        if (collision != empty) {
            fVersor versor = fVersor::reflector(collision);
            this->translator = versor(this->translator);
            if (this->bounce-- <= 0) {
                return false;
            }
        }
        
        return true;
    }
    
    
    void Projectile::buildVertices(std::vector<ObjectVertex> &vertices) const {
        GLfloat invLatitude, invLongitude, stepPhi, stepTheta, radius, angle;
        std::vector<ObjectVertex> data;
        c3ga::Mvec<GLfloat> dual;
        glm::vec3 center;
        
        invLatitude = 1.f / this->latitude;
        invLongitude = 1.f / this->longitude;
        stepPhi = 2 * M_PIf32 * invLatitude;
        stepTheta = M_PIf32 * invLongitude;
        
        dual = !this->multivector;
        dual /= dual[c3ga::E0];
        radius = std::abs(std::sqrt(dual | dual));
        center = { dual[c3ga::E1], dual[c3ga::E2], dual[c3ga::E3] };
        for (GLsizei j = 0; j <= this->longitude; ++j) {
            angle = -M_PI_2f32 + j * stepTheta;
            GLfloat cosTheta = std::cos(angle);
            GLfloat sinTheta = std::sin(angle);
            
            for (GLsizei i = 0; i <= this->latitude; ++i) {
                ObjectVertex vertex {};
                
                vertex.texture.x = i * invLatitude;
                vertex.texture.y = 1.f - j * invLongitude;
                vertex.normal.x = std::sin(i * stepPhi) * cosTheta;
                vertex.normal.y = sinTheta;
                vertex.normal.z = std::cos(i * stepPhi) * cosTheta;
                vertex.position = radius * vertex.normal + center;
                
                data.push_back(vertex);
            }
        }
        
        for (GLuint j = 0; j < this->longitude; ++j) {
            GLuint offset = j * (this->latitude + 1);
            
            for (GLuint i = 0; i < this->latitude; ++i) {
                vertices.push_back(data[offset + i]);
                vertices.push_back(data[offset + (i + 1)]);
                vertices.push_back(data[offset + this->latitude + 1 + (i + 1)]);
                vertices.push_back(data[offset + i]);
                vertices.push_back(data[offset + this->latitude + 1 + (i + 1)]);
                vertices.push_back(data[offset + i + this->latitude + 1]);
            }
        }
    }
    
    
    void Projectile::fillGlData(const std::vector<ObjectVertex> &vertices) const {
        // Fill the VBO
        glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(ObjectVertex) * this->size, vertices.data(), GL_STATIC_DRAW);
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
    
    
    bool Projectile::update() {
        std::vector<ObjectVertex> vertices;
        
        this->multivector = this->translator(this->multivector);
        if (!checkCollision()) {
            return false;
        }
        buildVertices(vertices);
        this->size = vertices.size();
        fillGlData(vertices);
        
        return true;
    }
    
    
    void Projectile::render() const {
        this->shader->bindTexture(*this->texture);
        glBindVertexArray(this->vao);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(this->size));
        glBindVertexArray(0);
        this->shader->unbindTexture();
    }
}
