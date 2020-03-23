#include <c3ga/Tools.hpp>

#include <glm/gtc/type_ptr.hpp>

#include <app/object/Triangle.hpp>
#include <app/object/Arena.hpp>
#include <app/object/Square.hpp>
#include <app/Engine.hpp>


namespace app::object {
    
    void Arena::xRing(const std::shared_ptr<shader::ShaderTexture> &shader,
                      const std::shared_ptr<shader::Texture> &texture) {
        fVersor rotorX90 = fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(90.));
        fVersor motor = (
                fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(45.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        
        c3ga::Mvec<GLfloat> a = motor(c3ga::point<GLfloat>(-0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> b = motor(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> c = motor(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        c3ga::Mvec<GLfloat> d = motor(c3ga::point<GLfloat>(-0.5f, 0.5f, 0));
        auto backTop = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("backTop", backTop);
        
        a = rotorX90(a);
        b = rotorX90(b);
        c = rotorX90(c);
        d = rotorX90(d);
        auto frontTop = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("frontTop", frontTop);
        
        a = rotorX90(a);
        b = rotorX90(b);
        c = rotorX90(c);
        d = rotorX90(d);
        auto frontBottom = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("frontBottom", frontBottom);
        
        a = rotorX90(a);
        b = rotorX90(b);
        c = rotorX90(c);
        d = rotorX90(d);
        auto backBottom = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("backBottom", backBottom);
    }
    
    
    void Arena::yRing(const std::shared_ptr<shader::ShaderTexture> &shader,
                      const std::shared_ptr<shader::Texture> &texture) {
        fVersor rotorY90 = fVersor::rotor(c3ga::e13<GLfloat>(), glm::radians(90.));
        fVersor motor = (
                fVersor::rotor(c3ga::e13<GLfloat>(), glm::radians(45.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        
        c3ga::Mvec<GLfloat> a = motor(c3ga::point<GLfloat>(-0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> b = motor(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> c = motor(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        c3ga::Mvec<GLfloat> d = motor(c3ga::point<GLfloat>(-0.5f, 0.5f, 0));
        auto backRight = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("backRight", backRight);
        
        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        d = rotorY90(d);
        auto frontRight = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("frontRight", frontRight);
        
        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        d = rotorY90(d);
        auto frontLeft = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("frontLeft", frontLeft);
        
        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        d = rotorY90(d);
        auto backLeft = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("backLeft", backLeft);
    }
    
    
    void Arena::zRing(const std::shared_ptr<shader::ShaderTexture> &shader,
                      const std::shared_ptr<shader::Texture> &texture) {
        fVersor rotorZ90 = fVersor::rotor(c3ga::e12<GLfloat>(), glm::radians(90.));
        fVersor motor = (
                fVersor::rotor(c3ga::e12<GLfloat>(), glm::radians(45.))
                * fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(90.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        
        c3ga::Mvec<GLfloat> a = motor(c3ga::point<GLfloat>(-0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> b = motor(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> c = motor(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        c3ga::Mvec<GLfloat> d = motor(c3ga::point<GLfloat>(-0.5f, 0.5f, 0));
        auto rightTop = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("rightTop", rightTop);
        
        a = rotorZ90(a);
        b = rotorZ90(b);
        c = rotorZ90(c);
        d = rotorZ90(d);
        auto rightBottom = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("rightBottom", rightBottom);

        a = rotorZ90(a);
        b = rotorZ90(b);
        c = rotorZ90(c);
        d = rotorZ90(d);
        auto leftBottom = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("leftBottom", leftBottom);

        a = rotorZ90(a);
        b = rotorZ90(b);
        c = rotorZ90(c);
        d = rotorZ90(d);
        auto leftTop = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("leftTop", leftTop);
    }
    
    
    void Arena::fillRings(const std::shared_ptr<shader::ShaderTexture> &shader,
                          const std::shared_ptr<shader::Texture> &texture) {
        fVersor rotorX90 = fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(90.));
        fVersor rotorZ90 = fVersor::rotor(c3ga::e12<GLfloat>(), glm::radians(90.));
        fVersor rotorZ180 = fVersor::rotor(c3ga::e12<GLfloat>(), glm::radians(180.));
        fVersor motor = fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>());
        
        c3ga::Mvec<GLfloat> a = motor(c3ga::point<GLfloat>(-0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> b = motor(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> c = motor(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        c3ga::Mvec<GLfloat> d = motor(c3ga::point<GLfloat>(-0.5f, 0.5f, 0));
        auto back = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("back", back);

        a = rotorX90(a);
        b = rotorX90(b);
        c = rotorX90(c);
        d = rotorX90(d);
        auto top = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("top", top);

        a = rotorX90(a);
        b = rotorX90(b);
        c = rotorX90(c);
        d = rotorX90(d);
        auto front = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("front", front);

        a = rotorX90(a);
        b = rotorX90(b);
        c = rotorX90(c);
        d = rotorX90(d);
        auto bottom = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("bottom", bottom);

        a = rotorZ90(a);
        b = rotorZ90(b);
        c = rotorZ90(c);
        d = rotorZ90(d);
        auto right = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("right", right);

        a = rotorZ180(a);
        b = rotorZ180(b);
        c = rotorZ180(c);
        d = rotorZ180(d);
        auto left = std::make_shared<object::Square>(a, b, c, d, shader, texture);
        this->rendered.emplace("left", left);
    }
    
    
    void Arena::triangles(const std::shared_ptr<shader::ShaderTexture> &shader,
                          const std::shared_ptr<shader::Texture> &black,
                          const std::shared_ptr<shader::Texture> &white) {
        fVersor rotorY90 = fVersor::rotor(c3ga::e13<GLfloat>(), glm::radians(90.));
        fVersor rotorX90 = fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(90.));
        
        fVersor motor1 = (
                fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(45.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        fVersor motor2 = (
                fVersor::rotor(c3ga::e13<GLfloat>(), glm::radians(45.))
                * fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(45.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        
        c3ga::Mvec<GLfloat> a = motor1(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> b = motor2(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        c3ga::Mvec<GLfloat> c = motor1(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        auto backRightTop = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("backRightTop", backRightTop);
        
        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        auto frontRightTop = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("frontRightTop", frontRightTop);
        
        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        auto frontLeftTop = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("frontLeftTop", frontLeftTop);
        
        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        auto backLeftTop = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("backLeftTop", backLeftTop);
        
        motor1 = (
                fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(-45.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        motor2 = (
                fVersor::rotor(c3ga::e13<GLfloat>(), glm::radians(45.))
                * fVersor::rotor(c3ga::e23<GLfloat>(), glm::radians(-45.))
                * fVersor::translator(-(SQRT_1_2 + 0.5f) * c3ga::e3<GLfloat>())
        );
        a = motor1(c3ga::point<GLfloat>(0.5f, -0.5f, 0));
        b = motor2(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        c = motor1(c3ga::point<GLfloat>(0.5f, 0.5f, 0));
        auto backRightBottom = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("backRightBottom", backRightBottom);

        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        auto frontRightBottom = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("frontRightBottom", frontRightBottom);

        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        auto frontLeftBottom = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("frontLeftBottom", frontLeftBottom);

        a = rotorY90(a);
        b = rotorY90(b);
        c = rotorY90(c);
        auto backLeftBottom = std::make_shared<object::Triangle>(a, b, c, shader, white);
        this->rendered.emplace("backLeftBottom", backLeftBottom);
    }
    
    
    Arena::Arena(std::shared_ptr<shader::ShaderTexture> shader) :
            shader(std::move(shader)), modified(true) {
        auto black = std::make_shared<shader::Texture>(misc::Image::loadPNG("../assets/black.png"));
        auto white = std::make_shared<shader::Texture>(misc::Image::loadPNG("../assets/white.png"));
        
        this->xRing(this->shader, black);
        this->yRing(this->shader, black);
        this->zRing(this->shader, black);
        this->fillRings(this->shader, white);
        this->triangles(this->shader, black, white);
    }
    
    
    void object::Arena::transform(const Versor<GLfloat> &versor) {
        std::for_each(this->rendered.begin(), this->rendered.end(), [&versor](auto r) { r.second->transform(versor); });
        this->modified = true;
    }
    
    
    c3ga::Mvec<GLfloat> Arena::collide(c3ga::Mvec<GLfloat> mvec) {
        c3ga::Mvec<GLfloat> empty = c3ga::Mvec<GLfloat>(), collision;
    
        for (const auto &r: this->rendered) {
            collision = r.second->collide(mvec);
            if (collision != empty) {
                return collision;
            }
        }
    
        return empty;
    }
    
    
    bool Arena::update() {
        if (!this->modified) {
            return true;
        }
        
        std::for_each(this->rendered.begin(), this->rendered.end(), [](auto r) { r.second->update(); });
        this->modified = false;
        return true;
    }
    
    
    void Arena::render() const {
        std::for_each(this->rendered.cbegin(), this->rendered.cend(), [](auto r) { r.second->render(); });
    }
}
