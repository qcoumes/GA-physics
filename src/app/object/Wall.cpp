#include <c3ga/Tools.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <app/object/Wall.hpp>
#include <app/object/ObjectVertex.hpp>
#include <app/Engine.hpp>


namespace app::object {
    
    Wall::Wall() :
            shader(shader::ShaderTexture("../shader/default.vs.glsl", "../shader/default.fs.glsl")),
            texture(shader::Texture(misc::Image::loadPNG("../assets/white4x4.png"))), modified(true) {
        
        this->a = c3ga::point<GLfloat>(0, 0, 0);
        this->b = c3ga::point<GLfloat>(1, 0, 0);
        this->c = c3ga::point<GLfloat>(1, 1, 0);
        this->d = c3ga::point<GLfloat>(0, 1, 0);
        this->plane = a ^ b ^ c ^ c3ga::ei<GLfloat>();
        
        glGenBuffers(1, &this->vbo);
        glGenVertexArrays(1, &this->vao);
        
        this->shader.addUniform("uMV", shader::UNIFORM_MATRIX_4F);
        this->shader.addUniform("uMVP", shader::UNIFORM_MATRIX_4F);
        this->shader.addUniform("uNormal", shader::UNIFORM_MATRIX_4F);
        this->shader.addUniform("uLightPosition", shader::UNIFORM_3_F);
        this->shader.addUniform("uLightColor", shader::UNIFORM_3_F);
        this->shader.addUniform("uLightDirIntensity", shader::UNIFORM_1_F);
        this->shader.addUniform("uLightAmbIntensity", shader::UNIFORM_1_F);
        GLfloat dirIntensity = 1.f, ambIntensity = 0.4f;
        this->shader.use();
        this->shader.loadUniform("uLightPosition", glm::value_ptr(glm::vec3(100, 100, -100)));
        this->shader.loadUniform("uLightColor", glm::value_ptr(glm::vec3(1, 1, 1)));
        this->shader.loadUniform("uLightDirIntensity", &dirIntensity);
        this->shader.loadUniform("uLightAmbIntensity", &ambIntensity);
        this->shader.stop();
    }
    
    
    Wall::~Wall() {
        glDeleteBuffers(1, &this->vbo);
        glDeleteVertexArrays(1, &this->vao);
    }
    
    
    void Wall::transform(const Versor<GLfloat> &versor) {
        std::cout << this->c[c3ga::E1] << " " << this->c[c3ga::E2] << " " << this->c[c3ga::E3] << std::endl;
        this->a = versor(this->a);
        this->b = versor(this->b);
        this->c = versor(this->c);
        this->d = versor(this->d);
        this->plane = a ^ b ^ c ^ c3ga::ei<GLfloat>();
        
        std::cout << this->c[c3ga::E1] << " " << this->c[c3ga::E2] << " " << this->c[c3ga::E3] << std::endl;
        
        modified = true;
    }
    
    
    void Wall::update() {
        glm::vec3 avec = { this->a[c3ga::E1], this->a[c3ga::E2], this->a[c3ga::E3] };
        glm::vec3 bvec = { this->b[c3ga::E1], this->b[c3ga::E2], this->b[c3ga::E3] };
        glm::vec3 cvec = { this->c[c3ga::E1], this->c[c3ga::E2], this->c[c3ga::E3] };
        glm::vec3 dvec = { this->d[c3ga::E1], this->d[c3ga::E2], this->d[c3ga::E3] };
        
        ObjectVertex vertices[6] = {
                ObjectVertex(avec, { 0, 0, 1 }, { 1, 1 }),
                ObjectVertex(bvec, { 0, 0, 1 }, { 0, 1 }),
                ObjectVertex(cvec, { 0, 0, 1 }, { 0, 0 }),
                ObjectVertex(cvec, { 0, 0, 1 }, { 0, 0 }),
                ObjectVertex(dvec, { 0, 0, 1 }, { 1, 0 }),
                ObjectVertex(avec, { 0, 0, 1 }, { 1, 1 }),
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
    }
    
    
    void Wall::render() const {
        Engine *engine = Engine::getInstance();
        glm::mat4 MVMatrix = engine->camera->getViewMatrix();
        glm::mat4 MVPMatrix = engine->camera->getProjMatrix() * MVMatrix;
        glm::mat4 normalMatrix = glm::transpose(glm::inverse(MVMatrix));
        
        this->shader.use();
        this->shader.loadUniform("uMV", glm::value_ptr(MVMatrix));
        this->shader.loadUniform("uMVP", glm::value_ptr(MVPMatrix));
        this->shader.loadUniform("uNormal", glm::value_ptr(normalMatrix));
        this->shader.bindTexture(this->texture);
        glBindVertexArray(this->vao);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(6));
        glBindVertexArray(0);
        this->shader.unbindTexture();
        this->shader.stop();
    }
}
