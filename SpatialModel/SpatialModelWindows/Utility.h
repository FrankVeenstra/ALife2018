#pragma once

#include <SFML/Graphics.hpp>

using namespace std ;
using namespace sf ;

float clamp(float value, float min, float max) ;
bool mouseClickedOnce(Mouse::Button const& button, bool &clickedLastFrame) ;
