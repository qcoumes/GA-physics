#include <c3ga/Tools.hpp>
#include <utility>
#include <glm/vec3.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <app/object/Square.hpp>
#include <app/object/ObjectVertex.hpp>


namespace app::object {
    
    Square::Square(c3ga::Mvec<GLfloat> a, c3ga::Mvec<GLfloat> b, c3ga::Mvec<GLfloat> c, c3ga::Mvec<GLfloat> d,
                   std::shared_ptr<shader::ShaderTexture> shader, std::shared_ptr<shader::Texture> texture) :
            shader(std::move(shader)), texture(std::move(texture)), a(std::move(a)), b(std::move(b)), c(std::move(c)),
            d(std::move(d)), modified(true) {
        
        this->plane = this->a ^ this->b ^ this->c ^ c3ga::ei<GLfloat>();
        this->normal = plane * this->a;
        
        glGenBuffers(1, &this->vbo);
        glGenVertexArrays(1, &this->vao);
    }
    
    
    Square::~Square() {
        glDeleteBuffers(1, &this->vbo);
        glDeleteVertexArrays(1, &this->vao);
    }
    
    
    void Square::transform(const Versor<GLfloat> &versor) {
        this->a = versor(this->a);
        this->b = versor(this->b);
        this->c = versor(this->c);
        this->d = versor(this->d);
        this->plane = this->a ^ this->b ^ this->c ^ c3ga::ei<GLfloat>();
        this->normal = this->plane * this->a;
    
        this->modified = true;
    }
    
    
    c3ga::Mvec<GLfloat> Square::collide(c3ga::Mvec<GLfloat> mvec) {
        c3ga::Mvec<GLfloat> collision = !mvec  | this->plane;
        if ((collision | collision) >= 0.f) {
            return this->plane;
        }
        return c3ga::Mvec<GLfloat>();
    }
    
    
    bool Square::update() {
        if (!this->modified) {
            return true;
        }
        
        c3ga::Mvec<GLfloat> orientation = (c3ga::ei<GLfloat>() | (!this->normal)) ^ c3ga::ei<GLfloat>();
        orientation = orientation | c3ga::e0<GLfloat>();
        orientation /= orientation.norm();
        glm::vec3 normal = { orientation[c3ga::E1], orientation[c3ga::E2], orientation[c3ga::E3] };
        
        c3ga::Mvec<GLfloat> a0 = this->a / this->a[c3ga::E0];
        c3ga::Mvec<GLfloat> b0 = this->b / this->b[c3ga::E0];
        c3ga::Mvec<GLfloat> c0 = this->c / this->c[c3ga::E0];
        c3ga::Mvec<GLfloat> d0 = this->d / this->d[c3ga::E0];
        glm::vec3 aPos = { a0[c3ga::E1], a0[c3ga::E2], a0[c3ga::E3] };
        glm::vec3 bPos = { b0[c3ga::E1], b0[c3ga::E2], b0[c3ga::E3] };
        glm::vec3 cPos = { c0[c3ga::E1], c0[c3ga::E2], c0[c3ga::E3] };
        glm::vec3 dPos = { d0[c3ga::E1], d0[c3ga::E2], d0[c3ga::E3] };
        
        ObjectVertex vertices[6] = {
                ObjectVertex(aPos, normal, { 1, 1 }),
                ObjectVertex(bPos, normal, { 0, 1 }),
                ObjectVertex(cPos, normal, { 0, 0 }),
                ObjectVertex(cPos, normal, { 0, 0 }),
                ObjectVertex(dPos, normal, { 1, 0 }),
                ObjectVertex(aPos, normal, { 1, 1 }),
        };
        
        // Fill the VBO
        glBindBuffer(GL_ARRAY_BUFFER, this->vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
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
        
        this->modified = false;
        return true;
    }
    
    
    void Square::render() const {
        this->shader->bindTexture(*this->texture);
        glBindVertexArray(this->vao);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(6));
        glBindVertexArray(0);
        this->shader->unbindTexture();
    }
}
