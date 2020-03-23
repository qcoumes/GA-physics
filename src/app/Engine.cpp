#include <chrono>
#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>

#include <imgui/imgui.h>
#include <glm/gtc/type_ptr.hpp>

#include <tool/ImGuiHandler.hpp>
#include <app/Engine.hpp>
#include <app/Config.hpp>
#include <app/Stats.hpp>
#include <app/object/Arena.hpp>
#include <app/object/Projectile.hpp>


namespace app {
    
    Engine *Engine::getInstance() {
        static Engine engine;
        return &engine;
    }
    
    
    void Engine::init() {
        this->window = std::make_unique<tool::Window>("OpenGL");
        
        GLenum glewInitError = glewInit();
        if (GLEW_OK != glewInitError) {
            throw std::runtime_error(reinterpret_cast<const char *>(glewGetErrorString(glewInitError)));
        }
        
        this->lastTick = std::chrono::steady_clock::now();
        this->camera = std::make_unique<tool::Camera>();
        //        this->camera->moveForward(-5);
        this->input = std::make_unique<tool::Input>();
        this->imGui = std::make_unique<tool::ImGuiHandler>(this->window->getWindow(), this->window->getContext());
        this->projTexture = std::make_shared<shader::Texture>(misc::Image::loadPNG("../assets/tennis.png"));
        
        Config::getInstance()->init(*this->window, *this->camera);
        
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        this->shader = std::make_shared<shader::ShaderTexture>(
                "../shader/default.vs.glsl", "../shader/default.fs.glsl"
        );
        this->shader->addUniform("uMV", shader::UNIFORM_MATRIX_4F);
        this->shader->addUniform("uMVP", shader::UNIFORM_MATRIX_4F);
        this->shader->addUniform("uNormal", shader::UNIFORM_MATRIX_4F);
        this->shader->addUniform("uLightPosition", shader::UNIFORM_3_F);
        this->shader->addUniform("uLightColor", shader::UNIFORM_3_F);
        this->shader->addUniform("uLightDirIntensity", shader::UNIFORM_1_F);
        this->shader->addUniform("uLightAmbIntensity", shader::UNIFORM_1_F);
        GLfloat dirIntensity = 1.f, ambIntensity = 0.4f;
        this->shader->use();
        this->shader->loadUniform("uLightPosition", glm::value_ptr(glm::vec3(100, 100, -100)));
        this->shader->loadUniform("uLightColor", glm::value_ptr(glm::vec3(1, 1, 1)));
        this->shader->loadUniform("uLightDirIntensity", &dirIntensity);
        this->shader->loadUniform("uLightAmbIntensity", &ambIntensity);
        this->shader->stop();
        
        this->rendered.emplace("Arena", std::make_shared<object::Arena>(shader));
        this->rendered["Arena"]->transform(fVersor::dilator(30.f));
    }
    
    
    GLboolean Engine::tick() {
        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - this->lastTick).count();
        
        if (duration >= MS_PER_TICK) {
            this->lastTick = now;
            this->tickCount++;
            this->tickSecond = (this->tickSecond + 1) % TICK_PER_SEC;
            return true;
        }
        
        return false;
    }
    
    
    c3ga::Mvec<GLfloat> Engine::collide(c3ga::Mvec<GLfloat> mvec) {
        c3ga::Mvec<GLfloat> empty = c3ga::Mvec<GLfloat>(), collision;
        
        for (const auto &r: this->rendered) {
            collision = r.second->collide(mvec);
            if (collision != empty) {
                return collision;
            }
        }
        
        return empty;
    }
    
    
    void Engine::update() {
        Config *config = Config::getInstance();
        GLfloat speed = config->getSpeed();
        SDL_Event event;
        this->input->reset();
        
        while (SDL_PollEvent(&event)) {
            switch (event.type) {
                case SDL_MOUSEMOTION:
                case SDL_MOUSEBUTTONDOWN:
                case SDL_MOUSEBUTTONUP:
                case SDL_MOUSEWHEEL:
                    if (this->imGui->wantCaptureMouse()) {
                        this->imGui->handleEvent(event);
                    }
                    else {
                        this->input->handleInput(event);
                    }
                    break;
                
                case SDL_KEYDOWN:
                case SDL_KEYUP:
                    if (this->imGui->wantCaptureKeyboard()) {
                        this->imGui->handleEvent(event);
                    }
                    else {
                        this->input->handleInput(event);
                    }
                    break;
                
                default:
                    this->input->handleInput(event);
            }
        }
        
        // Close app
        if (this->input->ended() || this->input->isReleasedKey(SDL_SCANCODE_ESCAPE)) {
            this->running = false;
        }
        
        // Movements
        if (this->input->isHeldKey(SDL_SCANCODE_A)) {
            this->camera->moveLeft(speed);
        }
        if (this->input->isHeldKey(SDL_SCANCODE_D)) {
            this->camera->moveLeft(-speed);
        }
        if (this->input->isHeldKey(SDL_SCANCODE_W)) {
            this->camera->moveForward(speed);
        }
        if (this->input->isHeldKey(SDL_SCANCODE_S)) {
            this->camera->moveForward(-speed);
        }
        if (this->input->isHeldKey(SDL_SCANCODE_SPACE)) {
            this->camera->moveUp(speed);
        }
        if (this->input->isHeldKey(SDL_SCANCODE_LCTRL)) {
            this->camera->moveUp(-speed);
        }
        
        // Camera rotation
        if (!config->getFreeMouse()) {
            
            /* Workaround to fix a bug where the mouse is able to leave the window even when capture when imgui is
             * enabled
             */
            config->setFreeMouse(true);
            config->setFreeMouse(false);
            ////////////////////////////
            
            glm::vec2 motion = this->input->getRelativeMotion();
            GLfloat sensitivity = config->getMouseSensitivity();
            this->camera->rotateLeft(-motion.x * sensitivity);
            this->camera->rotateUp(-motion.y * sensitivity);
            
            
            // Throw ball
            if (this->input->isPressedButton(SDL_BUTTON_LEFT)) {
                this->projectile.insert(std::make_shared<object::Projectile>(this->shader, this->projTexture));
            }
        }
        
        // Change projectile settings
        if (this->input->isHeldButton(SDL_BUTTON_RIGHT)) {
            config->setProjSize(config->getProjSize() + (this->input->getWheelMotion().y) * 0.5f);
        }
        else if (this->input->isHeldButton(SDL_BUTTON_MIDDLE)) {
            config->setProjBounce(config->getProjBounce() + static_cast<GLuint>(this->input->getWheelMotion().y));
        }
        else {
            config->setProjSpeed(config->getProjSpeed() + (this->input->getWheelMotion().y * 0.2f));
        }
        
        // Toggle free mouse
        if (this->input->isReleasedKey(SDL_SCANCODE_LALT)) {
            config->switchFreeMouse();
        }
        // Toggle debug
        if (this->input->isReleasedKey(SDL_SCANCODE_F10)) {
            config->switchDebug();
        }
        
        std::vector<std::shared_ptr<object::ObjectC3GA>> toDelete;
        std::for_each(this->rendered.begin(), this->rendered.end(), [](auto r) { r.second->update(); });
        std::for_each(
                this->projectile.begin(), this->projectile.end(),
                [&toDelete](auto p) {
                    if (!p->update()) {
                        toDelete.push_back(p);
                    }
                }
        );
        std::for_each(toDelete.begin(), toDelete.end(), [this](auto p) { this->projectile.erase(p); });
    }
    
    
    void Engine::debug() const {
        Config *config = Config::getInstance();
        Engine *engine = Engine::getInstance();
        glm::dvec3 position = engine->camera->getPosition();
        glm::dvec3 lookingAt = engine->camera->getFrontVector();
        std::stringstream ss;
        GLint length;
        
        this->imGui->newFrame();
        ImGui::SetNextWindowSize({ 600, 500 }, ImGuiCond_Once);
        
        ImGui::Begin("Debug");
        ImGui::Text("FPS: %d /", Stats::getInstance()->fps);
        ImGui::SameLine();
        ImGui::PushItemWidth(105);
        if (ImGui::BeginCombo("", config->getFramerateString().c_str())) {
            for (GLint i = 0; i <= Framerate::FRAMERATE_LAST; i++) {
                bool is_selected = (config->getFramerateOpt() == i);
                if (ImGui::Selectable(config->getFramerateString(i).c_str(), is_selected)) {
                    config->setFramerate(static_cast<Framerate>(i));
                }
                if (is_selected) {
                    ImGui::SetItemDefaultFocus();
                }
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        tool::ImGuiHandler::HelpMarker("May be overriden by your driver / windowing system");
        length = static_cast<GLint>(std::to_string(Config::TICK_PER_SEC).size());
        ss.str(std::string());
        ss << "Tick: " << std::setfill('0') << std::setw(length) << engine->tickSecond << "/" << Config::TICK_PER_SEC;
        ImGui::Text("%s", ss.str().c_str());
        ImGui::Text("Position: (%.2f, %.2f, %.2f)", position.x, position.y, position.z);
        ImGui::Text("Looking at: (%.2f, %.2f, %.2f)", lookingAt.x, lookingAt.y, lookingAt.z);
        ImGui::Dummy({ 0.0f, 6.0f });
        
        ImGui::Text("CONTROLES:");
        ImGui::BulletText("W, A, S, D, CTL, SPACE : Move the camera.");
        ImGui::BulletText("LEFT CLICK: Throw a ball.");
        ImGui::BulletText("MOUSE WHEEL: Change the speed of the ball.");
        ImGui::BulletText("RIGHT CLICK + MOUSE WHEEL: Change the size of the ball.");
        ImGui::BulletText("MIDDLE CLICK + MOUSE WHEEL: Change the number of bounce before a ball disapear.");
        ImGui::BulletText("LEFT ALT: Free / lock the mouse.");
        ImGui::Dummy({ 0.0f, 3.0f });
        ImGui::Text(
                "Once the mouse is freed, you can modify speed, size and bounce\n"
                "directly with the sliders below."
        );
        ImGui::Dummy({ 0.0f, 6.0f });
        
        if (ImGui::CollapsingHeader("Projectiles", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("Changes are only applied to new balls.");
            
            GLfloat speed = config->getProjSpeed(), size = config->getProjSize();
            GLint bounce = config->getProjBounce();
            if (ImGui::SliderFloat("Speed", &speed, 0.3f, 3.f)) {
                config->setProjSpeed(speed);
            }
            if (ImGui::SliderFloat("Size", &size, 0.5f, 10.f)) {
                config->setProjSize(size);
            }
            if (ImGui::SliderInt("Bounce", &bounce, 1, 20)) {
                config->setProjBounce(bounce);
            }
            ImGui::Text(
                    "Note that if a ball is too small and too fast, it could move through\n"
                    "a wall before the collision is checked."
            );
        }
        
        if (ImGui::CollapsingHeader("Hardware & Driver")) {
            GLint offset = 120;
            ImGui::Text("CPU:");
            ImGui::SameLine(offset);
            ImGui::Text("%s", config->getCpuInfo().c_str());
            ImGui::Separator();
            
            ImGui::Text("GPU:");
            ImGui::SameLine(offset);
            ImGui::Text("%s", config->getGPUInfo().c_str());
            ImGui::Text("Driver:");
            ImGui::SameLine(offset);
            ImGui::Text("%s", config->getGPUDriver().c_str());
            if (ImGui::TreeNode("Extensions")) {
                for (const std::string &s : config->getGPUExtensions()) {
                    ImGui::BulletText("%s", s.c_str());
                }
                ImGui::TreePop();
            }
            ImGui::Separator();
            
            ImGui::Text("GLEW Version:");
            ImGui::SameLine(offset);
            ImGui::Text("%s", config->getGlewVersion().c_str());
        }
        
        if (ImGui::CollapsingHeader("Settings")) {
            bool faceCulling = config->getFaceCulling();
            ImGui::Checkbox("Face Culling", &faceCulling);
            if (faceCulling != config->getFaceCulling()) {
                config->switchFaceCulling();
            }
        }
        ImGui::End();
        this->imGui->render();
    }
    
    
    void Engine::_render() const {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        Engine *engine = Engine::getInstance();
        glm::mat4 MVMatrix = engine->camera->getViewMatrix();
        glm::mat4 MVPMatrix = engine->camera->getProjMatrix() * MVMatrix;
        glm::mat4 normalMatrix = glm::transpose(glm::inverse(MVMatrix));
        
        this->shader->use();
        this->shader->loadUniform("uMV", glm::value_ptr(MVMatrix));
        this->shader->loadUniform("uMVP", glm::value_ptr(MVPMatrix));
        this->shader->loadUniform("uNormal", glm::value_ptr(normalMatrix));
        std::for_each(this->rendered.cbegin(), this->rendered.cend(), [](const auto &r) { r.second->render(); });
        std::for_each(this->projectile.cbegin(), this->projectile.cend(), [](const auto &p) { p->render(); });
        this->shader->stop();
        
        if (Config::getInstance()->getDebug()) {
            this->debug();
        }
        
        this->window->refresh();
    }
    
    
    void Engine::render() const {
        static std::chrono::steady_clock::time_point cmptStart = std::chrono::steady_clock::now();
        static std::chrono::steady_clock::time_point limiterStart = cmptStart;
        static GLuint fps = 0;
        
        std::chrono::steady_clock::time_point now;
        double duration;
        
        if (Config::getInstance()->getFramerate() > 0) { // Capped framerate
            now = std::chrono::steady_clock::now();
            duration = std::chrono::duration_cast<std::chrono::microseconds>(now - limiterStart).count();
            
            if (duration >= Config::getInstance()->getFramerateInv()) {
                this->_render();
                limiterStart = now;
                fps++;
            }
        }
        else { // Uncapped framerate
            this->_render();
            fps++;
        }
        
        // Computing FPS
        now = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::seconds>(now - cmptStart).count();
        if (duration >= 1.) {
            Stats::getInstance()->fps = static_cast<GLuint>(fps);
            fps = 0;
            cmptStart = now;
        }
    }
    
    
    void Engine::cleanup() {
    }
    
    
    GLboolean Engine::isRunning() {
        return this->running;
    }
}
