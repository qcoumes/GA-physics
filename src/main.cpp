#include <app/Engine.hpp>


using namespace app;


int main(int argc, char **argv) {
    Engine *game = Engine::getInstance();
    
    std::cout << sizeof(c3ga::Mvec<GLfloat>) << std::endl;
    
    game->init();
    game->update();
    
    while (game->isRunning()) {
        if (game->tick()) {
            game->update();
        }
        
        game->render();
    }
    
    game->cleanup();
    
    return EXIT_SUCCESS;
}

