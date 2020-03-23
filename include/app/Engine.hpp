#ifndef OPENGL_ENGINE_HPP
#define OPENGL_ENGINE_HPP

#include <unordered_set>
#include <unordered_map>
#include <cstdint>
#include <chrono>
#include <memory>

#include <GL/glew.h>

#include <shader/ShaderTexture.hpp>
#include <shader/Texture.hpp>
#include <tool/Window.hpp>
#include <tool/Camera.hpp>
#include <tool/Input.hpp>
#include <tool/ImGuiHandler.hpp>
#include <app/object/ObjectC3GA.hpp>


namespace app {
    
    class Engine : public misc::ISingleton {
        private:
            static constexpr GLuint64 TICK_PER_SEC = 60;
            static constexpr GLdouble MS_PER_TICK = 1. / TICK_PER_SEC * 1000.;
            
            std::shared_ptr<shader::ShaderTexture> shader;
            std::shared_ptr<shader::Texture> projTexture;
            
            std::chrono::steady_clock::time_point lastTick;
            GLboolean running = true;
            GLuint tickSecond = 0;
            GLuint tickCount = 0;
        
        public:
            std::unordered_map<std::string, std::shared_ptr<object::ObjectC3GA>> rendered;
            std::unordered_set<std::shared_ptr<object::ObjectC3GA>> projectile;
            std::unique_ptr<tool::ImGuiHandler> imGui = nullptr;
            std::unique_ptr<tool::Window> window = nullptr;
            std::unique_ptr<tool::Camera> camera = nullptr;
            std::unique_ptr<tool::Input> input = nullptr;
        
        private:
            
            void debug() const;
            
            void _render() const;
            
            Engine() = default;
        
        public:
            
            static Engine *getInstance();
            
            void init();
            
            GLboolean tick();
        
            [[nodiscard]] c3ga::Mvec<GLfloat> collide(c3ga::Mvec<GLfloat> mvec);
            
            void update();
            
            void render() const;
            
            void cleanup();
            
            [[nodiscard]] GLboolean isRunning();
            
    };
}

#endif // OPENGL_ENGINE_HPP
